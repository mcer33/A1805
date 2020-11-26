import numpy as np
from numpy import linalg as LA
import math


def generateCoordinates(num, rad):
  coor = np.empty([num, 2])
  deg_step = 2*math.pi/num
  for i in range(num):
    coor[i, :] = [math.cos(i*deg_step), math.sin(i*deg_step)]
  coor *= rad
  return coor

def initDipoles(coor, mag):
  n_coor = len(coor)
  dip = np.empty([n_coor, 2])
  for i in range(n_coor):
    dip[i,:] = -coor[i] / LA.norm(coor[i])
  dip *= mag
  return dip

def createDipNeighbors(num):
  neighbors=[]
  for i in range(num):
    prev = num-1 if i == 0 else i-1
    neighbors.append([prev, i])
    next = 0 if i == num-1 else i+1
    neighbors.append([next, i])
  return neighbors

def createChainNeighbors(num):
  neighbors = []
  for i in range(num):
    next = 0 if i == num-1 else i+1
    neighbors.append([next, i])
  return neighbors

def applyFieldForce():
  i=0
  for dip in dipoles:
    tForces[i] += -Field*dip[0] #corresponds to -Y electric field direction
    i+=1

def getAngle(vec1,vec2):
  det = vec1[0] * vec2[1] - vec1[1] * vec2[0]
  dot = vec1[0] * vec2[0] + vec1[1] * vec2[1]
  angle = -math.atan2(det, dot)
  return angle

def clearForces(forceArr):
  forceArr.fill(0)

def applyDip2DipTorque(neigh):
  for pair in neigh:
    indx1 = pair[0]
    indx2 = pair[1]
    coor1 = coordinates[indx1]
    coor2 = coordinates[indx2]
    d_vec = coor2-coor1
    r = LA.norm(d_vec)
    phi1 = getAngle(dipoles[indx1],d_vec)
    phi2 = getAngle(dipoles[indx2],d_vec)
    tForces[indx2] += -dip_const*dip_mag**2/r**3*(-math.cos(phi1)*math.sin(phi2)-1/2*math.sin(phi1)*math.cos(phi2))

def generateJointPosition(pos, dir):
  sdir = np.array([-dir[1],dir[0]])
  sdir = sdir/LA.norm(sdir)
  donor = pos+p_radius*sdir
  acceptor = pos-p_radius*sdir
  return([donor,acceptor])

def applyCohesiveForce(neigh):
  for pair in neigh:
    indx1 = pair[0]
    indx2 = pair[1]
    coor1 = coordinates[indx1]
    coor2 = coordinates[indx2]
    dip1 = dipoles[indx1]
    dip2 = dipoles[indx2]
    joints_1 = generateJointPosition(coor1, dip1)
    joints_2 = generateJointPosition(coor2, dip2)
    r_sh = joints_2[1]-joints_1[0]
    mag_r = LA.norm(r_sh)
    force=-LJepsilon*((LJdelta/mag_r)**12 - (LJdelta/mag_r)**6)
    lForces[indx1] += force * r_sh/mag_r
    lForces[indx2] += -force * r_sh/mag_r

def rotateAntiClockwise(vec,ang_rad):
  x = math.cos(ang_rad)*vec[0] - math.sin(ang_rad)*vec[1]
  y = math.sin(ang_rad)*vec[0] + math.cos(ang_rad)*vec[1]
  return np.array([x,y])

def integrateStep(dStep):
  for i in range(num):
    beta = tForces[i]*dStep
    xy = rotateAntiClockwise(dipoles[i, :], beta)
    dipoles[i, :] = xy
    coordinates[i] += lForces[i] * d_step


####################################################################################
num = 13  # number of dipoles
coordinates = generateCoordinates(num, 90)
tForces = np.empty([num, 1])  # list of torque forces
lForces = np.empty([num, 2])  # list of linear forces
dip_mag = 15  # dipole moment magnitude (arbitrary units)
dip_const = 1 # multiplicative constant for dipole-dipole interaction
dipoles = initDipoles(coordinates, dip_mag)
p_radius = 18  # beads radius (arbitrary units)
LJepsilon = 15  # L-J Epsilon (arbitrary units)
LJdelta = 8  # L-J Delta (arbitrary units)
n_step = 10000  # total number of simulation steps
d_step = 0.01  # integration step (arbitrary units)
Field = 0.001  # field intensity (arbitrary units)
save_freq = 100  # snapshot output frequency
traj_path = "/home/jirka/WORK/DD/"  # path to save snapshots
dipNeighbors = createDipNeighbors(num)
chainNeighbors = createChainNeighbors(num)
startFieldTime = 1000  # which step to start apply the e-field force
####################################################################################

frame = 0
for step in range(n_step):
  if not step%save_freq:
    np.save(traj_path + "coor" + str(frame), coordinates)
    np.save(traj_path + "dip" + str(frame), dipoles)
    frame += 1
    print(step)

  clearForces(tForces)
  clearForces(lForces)
  if step > startFieldTime:
    applyFieldForce()
  applyDip2DipTorque(dipNeighbors)
  applyCohesiveForce(chainNeighbors)

  integrateStep(d_step)

  # print(step)

