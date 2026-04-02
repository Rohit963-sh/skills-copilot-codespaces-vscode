#uniaxial compressive strength (UCS) test
#gen cylindrical rock sample
#
#import the appropriate ESyS-Particle modules:
from esys.lsm import *
from esys.lsm.util import *
from esys.lsm.geometry import *
from WallLoader import WallLoaderRunnable
from ServoWallLoader import ServoWallLoaderRunnable
from outputData import *
from math import sqrt, pi
import numpy as np

#instantiate a simulation object and
#initialise the neighbour search algorithm:
sim = LsmMpi(numWorkerProcesses=8, mpiDimList=[2,2,2])
sim.initNeighbourSearch(
    particleType = "RotSphere",
    gridSpacing = 2.5*2,        #>=2.5*Rmax+0.2*Rmin
    verletDist = 0.1*2          #>=0.2*Rmin
)
#specify the spatial domain and direction of periodic boundaries:
domain = BoundingBox ( Vec3 (0,0,0), Vec3 (150,150,150) )
sim.setSpatialDomain (
    bBox = domain,
    circDimList = [False, False, False]
)
############################# constant parameters ##################################
rhos = 0.00098                          # particle density
beta = 0.08                             # damping coefficient

########################### geometry ##############################################
sampleHeight = 149.0                    # sample height (mm)
sampleArea   = 1.0                      # cross-sectional area (placeholder)

############################ physical properties ##################################
PoissonsRatio  = 0.33
YoungsModulus  = 15.0                   # Young's modulus (MPa or consistent units)
ShearModulus   = YoungsModulus / (2.0 * (1.0 + PoissonsRatio))

cohesion0      = 0.3                    # bond cohesion strength
tanAngle0      = 0.5                    # bond friction angle (tan)
bondTag        = 5                      # bond tag used in BrittleBeamPrms

frictionCoeff  = 0.1                    # particle-particle friction coefficient

kn = 15.0                               # normal stiffness of wall-ball contact

########################### simulation time control #################################
inter1000 = 1000
inter      = 1000
N          = 2000000                    # total simulation timesteps

######### UCS compression control ###############
strainrate0   = 0.03                    # wall velocity magnitude (mm/timestep)
startTime     = 10000
endTime       = N

##################################################
#specify the number of timesteps and the timestep increment:
sim.setNumTimeSteps(N)
sim.setTimeStepSize(0.006)

############################################################################
sim.readGeometry("ucs.geo")
Num = sim.getNumParticles()
print('Total number of particles in the model: ', Num)

############################################################################
#set particle density:
sim.setParticleDensity(
    tag=1,
    mask=-1,
    Density=rhos
)

#specify gravity (set to zero for UCS test):
sim.createInteractionGroup(
    GravityPrms(
        name="gravity",
        acceleration=Vec3(0.0, 0.0, 0.0)
    )
)

#add a bottom wall:
sim.createWall(
    name="bottom_wall",
    posn=Vec3(0.0, 0.0, 0.0),
    normal=Vec3(0.0, 1.0, 0.0)
)

#add a top wall:
sim.createWall(
    name="top_wall",
    posn=Vec3(0.0, sampleHeight, 0.0),
    normal=Vec3(0.0, -1.0, 0.0)
)
sim.createWall(
	name   = "left_wall",
	posn   = Vec3(20, 0.0, 0.0),
	normal = Vec3(1.0, 0.0, 0.0)
)
 
 
#add a right wall (expanded):
sim.createWall(
	name   = "right_wall",
	posn   = Vec3(150, 0.0, 0.0),
	normal = Vec3(-1.0, 0.0, 0.0)
)
#create rotational elastic-brittle bonds between particles:
pp_bonds = sim.createInteractionGroup(
    BrittleBeamPrms(
        name="bonded_grain",
        youngsModulus=YoungsModulus,
        poissonsRatio=PoissonsRatio,
        cohesion=cohesion0,
        tanAngle=tanAngle0,
        tag=bondTag
    )
)
fip=RotFrictionPrms( "friction", pi*YoungsModulus/2.0, frictionCoeff, frictionCoeff, pi*ShearModulus/2.0)    #friction

sim.createInteractionGroup(fip)
sim.createExclusion("bonded_grain","friction")


#add translational local damping:
sim.createInteractionGroup(
    LocalDampingPrms(
        name="damping1",
        viscosity=beta
    )
)

#add rotational local damping:
sim.createInteractionGroup(
    RotLocalDampingPrms(
        name="damping2",
        viscosity=beta
    )
)
sim.createInteractionGroup(
	NRotElasticWallPrms(
		name     = "left_repel",
		wallName = "left_wall",
		normalK  = 0.1
	)
)
 
 
#wall-particle contact stiffness — right wall:
sim.createInteractionGroup(
	NRotElasticWallPrms(
		name     = "right_repel",
		wallName = "right_wall",
		normalK  = 0.1
	)
)
# FIX 3: add front and back walls in Z — your geometry runs from Z=0 to Z≈130
sim.createWall(name="front_wall", posn=Vec3(0.0, 0.0, 0.0),   normal=Vec3(0.0, 0.0,  1.0))
sim.createWall(name="back_wall",  posn=Vec3(0.0, 0.0, 140.0), normal=Vec3(0.0, 0.0, -1.0))

sim.createInteractionGroup(NRotElasticWallPrms(name="front_repel", wallName="front_wall", normalK=0.1))
sim.createInteractionGroup(NRotElasticWallPrms(name="back_repel",  wallName="back_wall",  normalK=0.1))
#wall-particle contact stiffness — bottom wall:
sim.createInteractionGroup(
    NRotElasticWallPrms(
        name="bow_repel",
        wallName="bottom_wall",
        normalK=0.1
    )
)

#wall-particle contact stiffness — top wall:
sim.createInteractionGroup(
    NRotElasticWallPrms(
        name="top_repel",
        wallName="top_wall",
        normalK=0.1
    )
)

#add a wall loader to move the top wall downward (compression):
wall_loader = WallLoaderRunnable(
    sim,
    "top_wall",
    Vec3(0.0, -strainrate0, 0.0),       # downward velocity
    startTime,
    endTime
)
sim.addPostTimeStepRunnable(wall_loader)

#add a CheckPointer to store simulation snapshots:
sim.createCheckPointer(
    CheckPointPrms(
        fileNamePrefix="txt/UCS",
        beginTimeStep=0,
        endTimeStep=N,
        timeStepIncr=inter
    )
)

#create a FieldSaver for kinetic energy:
sim.createFieldSaver(
    ParticleScalarFieldSaverPrms(
        fieldName="e_kin",
        fileName="data1/ekin.dat",
        fileFormat="SUM",
        beginTimeStep=0,
        endTimeStep=N,
        timeStepIncr=inter1000
    )
)

#create a FieldSaver for bond count:
sim.createFieldSaver(
    InteractionScalarFieldSaverPrms(
        interactionName="bonded_grain",
        fieldName="count",
        fileName="data1/nbonds.dat",
        fileFormat="SUM",
        beginTimeStep=0,
        endTimeStep=N,
        timeStepIncr=inter1000
    )
)

#create a FieldSaver for wall positions:
posn_saver = WallVectorFieldSaverPrms(
    wallName=["top_wall", "bottom_wall"],
    fieldName="Position",
    fileName="data1/out_wallPosition.dat",
    fileFormat="RAW_SERIES",
    beginTimeStep=0,
    endTimeStep=N,
    timeStepIncr=inter1000
)
sim.createFieldSaver(posn_saver)

#create a FieldSaver for wall forces:
force_saver = WallVectorFieldSaverPrms(
    wallName=["top_wall", "bottom_wall"],
    fieldName="Force",
    fileName="data1/out_wallForce.dat",
    fileFormat="RAW_SERIES",
    beginTimeStep=0,
    endTimeStep=N,
    timeStepIncr=inter1000
)
sim.createFieldSaver(force_saver)

#create a FieldSaver for friction interaction forces:
RAW2_NormalForcesaver = InteractionVectorFieldSaverPrms(
    interactionName="friction",
    fieldName="force",
    fileName="data1/out_RAW_WITH_IDNormalForce",
    fileFormat="RAW_WITH_ID",
    beginTimeStep=0,
    endTimeStep=N,
    timeStepIncr=inter
)
sim.createFieldSaver(RAW2_NormalForcesaver)

#create a FieldSaver for bond interaction forces:
RAW2_bonds_NormalForcesaver = InteractionVectorFieldSaverPrms(
    interactionName="bonded_grain",
    fieldName="force",
    fileName="data1/out_bonds_RAW_WITH_IDNormalForce",
    fileFormat="RAW_WITH_ID",
    beginTimeStep=0,
    endTimeStep=N,
    timeStepIncr=inter
)
sim.createFieldSaver(RAW2_bonds_NormalForcesaver)

print("Starting UCS simulation...")
#execute the simulation:
sim.run()
