import os
import sys
from time import sleep
'''
resultspath is wherever your text file of 
streamline locations is 
'''
resultspath = "Dropbox/CIS/CIS410/custom_project/"

sys.path.append(os.path.abspath(resultspath))

from results import *
'''
Or whever you put the aneurysm data set
'''
OpenDatabase("localhost:/Users/caseychesshir/Dropbox/CIS/CIS410/custom_project/aneurysm_tutorial_data/aneurysm_data/aneurysm*.silo database", 0)

def initializeMesh():


	AddPlot("Subset", "Mesh", 1, 1)
	SubsetAtts = SubsetAttributes()
	SubsetAtts.colorType = SubsetAtts.ColorByMultipleColors  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
	SubsetAtts.colorTableName = "Default"
	SubsetAtts.invertColorTable = 0
	SubsetAtts.filledFlag = 1
	SubsetAtts.legendFlag = 1
	SubsetAtts.lineStyle = SubsetAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
	SubsetAtts.lineWidth = 0
	SubsetAtts.singleColor = (0, 0, 0, 255)
	SubsetAtts.SetMultiColor(0, (255, 0, 0, 255))
	SubsetAtts.subsetNames = ("Whole mesh (Mesh)")
	SubsetAtts.subsetType = SubsetAtts.Mesh  # Domain, Group, Material, EnumScalar, Mesh, Unknown
	SubsetAtts.opacity = 0.25098
	SubsetAtts.wireframe = 0
	SubsetAtts.drawInternal = 0
	SubsetAtts.smoothingLevel = 0
	SubsetAtts.pointSize = 0.05
	SubsetAtts.pointType = SubsetAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
	SubsetAtts.pointSizeVarEnabled = 0
	SubsetAtts.pointSizeVar = "default"
	SubsetAtts.pointSizePixels = 2
	SetPlotOptions(SubsetAtts)

def initializeStreamline():
	#good starting spot: 3.5 3.25 5.4
	#also starting spot: 3.51 3.41 5.51 
	AddPlot("Streamline", "velocity", 1, 1)

	StreamlineAtts = StreamlineAttributes()
	StreamlineAtts.sourceType = StreamlineAtts.SpecifiedCircle  # SpecifiedPoint, SpecifiedPointList, SpecifiedLine, SpecifiedCircle, SpecifiedPlane, SpecifiedSphere, SpecifiedBox, Selection
	StreamlineAtts.pointSource = (0, 0, 0)
	StreamlineAtts.lineStart = (0, 0, 0)
	StreamlineAtts.lineEnd = (1, 0, 0)
	StreamlineAtts.planeOrigin = (0, 0, 0)
	StreamlineAtts.planeNormal = (0, 0, 1)
	StreamlineAtts.planeUpAxis = (0, 1, 0)
	StreamlineAtts.radius = 1
	StreamlineAtts.sphereOrigin = (0, 0, 0)
	StreamlineAtts.boxExtents = (0, 1, 0, 1, 0, 1)
	StreamlineAtts.useWholeBox = 1
	StreamlineAtts.pointList = (0, 0, 0, 1, 0, 0, 0, 1, 0)
	StreamlineAtts.sampleDensity0 = 2
	StreamlineAtts.sampleDensity1 = 2
	StreamlineAtts.sampleDensity2 = 2
	StreamlineAtts.coloringMethod = StreamlineAtts.ColorByTime  # Solid, ColorBySpeed, ColorByVorticity, ColorByLength, ColorByTime, ColorBySeedPointID, ColorByVariable, ColorByCorrelationDistance, ColorByNumberDomainsVisited
	StreamlineAtts.colorTableName = "Default"
	StreamlineAtts.singleColor = (0, 0, 0, 255)
	StreamlineAtts.legendFlag = 1
	StreamlineAtts.lightingFlag = 1
	StreamlineAtts.integrationDirection = StreamlineAtts.Forward  # Forward, Backward, Both
	StreamlineAtts.maxSteps = 1000
	StreamlineAtts.terminateByDistance = 0
	StreamlineAtts.termDistance = 10
	StreamlineAtts.terminateByTime = 0
	StreamlineAtts.termTime = 10
	StreamlineAtts.maxStepLength = 0.1
	StreamlineAtts.limitMaximumTimestep = 0
	StreamlineAtts.maxTimeStep = 0.1
	StreamlineAtts.relTol = 0.0001
	StreamlineAtts.absTolSizeType = StreamlineAtts.FractionOfBBox  # Absolute, FractionOfBBox
	StreamlineAtts.absTolAbsolute = 1e-06
	StreamlineAtts.absTolBBox = 1e-06
	StreamlineAtts.fieldType = StreamlineAtts.Default  # Default, FlashField, M3DC12DField, M3DC13DField, Nek5000Field, NIMRODField
	StreamlineAtts.fieldConstant = 1
	StreamlineAtts.velocitySource = (0, 0, 0)
	StreamlineAtts.integrationType = StreamlineAtts.DormandPrince  # Euler, Leapfrog, DormandPrince, AdamsBashforth, RK4, M3DC12DIntegrator
	StreamlineAtts.parallelizationAlgorithmType = StreamlineAtts.VisItSelects  # LoadOnDemand, ParallelStaticDomains, MasterSlave, VisItSelects
	StreamlineAtts.maxProcessCount = 10
	StreamlineAtts.maxDomainCacheSize = 3
	StreamlineAtts.workGroupSize = 32
	StreamlineAtts.pathlines = 0
	StreamlineAtts.pathlinesOverrideStartingTimeFlag = 0
	StreamlineAtts.pathlinesOverrideStartingTime = 0
	StreamlineAtts.pathlinesPeriod = 0
	StreamlineAtts.pathlinesCMFE = StreamlineAtts.POS_CMFE  # CONN_CMFE, POS_CMFE
	StreamlineAtts.coordinateSystem = StreamlineAtts.AsIs  # AsIs, CylindricalToCartesian, CartesianToCylindrical
	StreamlineAtts.phiScalingFlag = 0
	StreamlineAtts.phiScaling = 1
	StreamlineAtts.coloringVariable = ""
	StreamlineAtts.legendMinFlag = 0
	StreamlineAtts.legendMaxFlag = 0
	StreamlineAtts.legendMin = 0
	StreamlineAtts.legendMax = 1
	StreamlineAtts.displayBegin = 0
	StreamlineAtts.displayEnd = 1
	StreamlineAtts.displayBeginFlag = 0
	StreamlineAtts.displayEndFlag = 0
	StreamlineAtts.referenceTypeForDisplay = StreamlineAtts.Distance  # Distance, Time, Step
	StreamlineAtts.displayMethod = StreamlineAtts.Lines  # Lines, Tubes, Ribbons
	StreamlineAtts.tubeSizeType = StreamlineAtts.FractionOfBBox  # Absolute, FractionOfBBox
	StreamlineAtts.tubeRadiusAbsolute = 0.125
	StreamlineAtts.tubeRadiusBBox = 0.005
	StreamlineAtts.ribbonWidthSizeType = StreamlineAtts.FractionOfBBox  # Absolute, FractionOfBBox
	StreamlineAtts.ribbonWidthAbsolute = 0.125
	StreamlineAtts.ribbonWidthBBox = 0.01
	StreamlineAtts.lineWidth = 2
	StreamlineAtts.showSeeds = 1
	StreamlineAtts.seedRadiusSizeType = StreamlineAtts.FractionOfBBox  # Absolute, FractionOfBBox
	StreamlineAtts.seedRadiusAbsolute = 1
	StreamlineAtts.seedRadiusBBox = 0.015
	StreamlineAtts.showHeads = 0
	StreamlineAtts.headDisplayType = StreamlineAtts.Sphere  # Sphere, Cone
	StreamlineAtts.headRadiusSizeType = StreamlineAtts.FractionOfBBox  # Absolute, FractionOfBBox
	StreamlineAtts.headRadiusAbsolute = 0.25
	StreamlineAtts.headRadiusBBox = 0.02
	StreamlineAtts.headHeightRatio = 2
	StreamlineAtts.opacityType = StreamlineAtts.FullyOpaque  # FullyOpaque, Constant, Ramp, VariableRange
	StreamlineAtts.opacityVariable = ""
	StreamlineAtts.opacity = 1
	StreamlineAtts.opacityVarMin = 0
	StreamlineAtts.opacityVarMax = 1
	StreamlineAtts.opacityVarMinFlag = 0
	StreamlineAtts.opacityVarMaxFlag = 0
	StreamlineAtts.tubeDisplayDensity = 10
	StreamlineAtts.geomDisplayQuality = StreamlineAtts.Medium  # Low, Medium, High, Super
	StreamlineAtts.sampleDistance0 = 10
	StreamlineAtts.sampleDistance1 = 10
	StreamlineAtts.sampleDistance2 = 10
	StreamlineAtts.fillInterior = 1
	StreamlineAtts.randomSamples = 0
	StreamlineAtts.randomSeed = 0
	StreamlineAtts.numberOfRandomSamples = 1
	StreamlineAtts.forceNodeCenteredData = 0
	StreamlineAtts.issueTerminationWarnings = 1
	StreamlineAtts.issueStiffnessWarnings = 1
	StreamlineAtts.issueCriticalPointsWarnings = 1
	StreamlineAtts.criticalPointThreshold = 0.001
	StreamlineAtts.varyTubeRadius = StreamlineAtts.None  # None, Scalar
	StreamlineAtts.varyTubeRadiusFactor = 10
	StreamlineAtts.varyTubeRadiusVariable = ""
	StreamlineAtts.correlationDistanceAngTol = 5
	StreamlineAtts.correlationDistanceMinDistAbsolute = 1
	StreamlineAtts.correlationDistanceMinDistBBox = 0.005
	StreamlineAtts.correlationDistanceMinDistType = StreamlineAtts.FractionOfBBox  # Absolute, FractionOfBBox
	StreamlineAtts.selection = ""
	SetPlotOptions(StreamlineAtts)
	# The AnimationPlay RPC is not supported in the VisIt module so it will not be logged.



def initializeClip1():
	AddOperator("ExternalSurface", 1)
	AddOperator("Clip", 1)
	ClipAtts = ClipAttributes()
	ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
	ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
	ClipAtts.plane1Status = 1
	ClipAtts.plane2Status = 0
	ClipAtts.plane3Status = 0
	ClipAtts.plane1Origin = (0, 0, 0)
	ClipAtts.plane2Origin = (0, 0, 0)
	ClipAtts.plane3Origin = (0, 0, 0)
	ClipAtts.plane1Normal = (1, 0, 0)
	ClipAtts.plane2Normal = (0, 1, 0)
	ClipAtts.plane3Normal = (0, 0, 1)
	ClipAtts.planeInverse = 1
	ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1  # None, Plane1, Plane2, Plane3
	ClipAtts.center = (3.4, 2.65, 5.5)
	ClipAtts.radius = 0.34
	ClipAtts.sphereInverse = 0
	SetOperatorOptions(ClipAtts, 1)

def initializeClip2():
	AddOperator("Clip", 2)
	ClipAtts = ClipAttributes()
	ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
	ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
	ClipAtts.plane1Status = 1
	ClipAtts.plane2Status = 0
	ClipAtts.plane3Status = 0
	ClipAtts.plane1Origin = (0, 0, 0)
	ClipAtts.plane2Origin = (0, 0, 0)
	ClipAtts.plane3Origin = (0, 0, 0)
	ClipAtts.plane1Normal = (1, 0, 0)
	ClipAtts.plane2Normal = (0, 1, 0)
	ClipAtts.plane3Normal = (0, 0, 1)
	ClipAtts.planeInverse = 1
	ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1  # None, Plane1, Plane2, Plane3
	ClipAtts.center = (3.4, 4.45, 5.8)
	ClipAtts.radius = 0.2
	ClipAtts.sphereInverse = 0
	SetOperatorOptions(ClipAtts, 2)


def initializeClip3():
	AddOperator("Clip", 3)
	ClipAtts = ClipAttributes()
	ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
	ClipAtts.funcType = ClipAtts.Sphere  # Plane, Sphere
	ClipAtts.plane1Status = 1
	ClipAtts.plane2Status = 0
	ClipAtts.plane3Status = 0
	ClipAtts.plane1Origin = (0, 0, 0)
	ClipAtts.plane2Origin = (0, 0, 0)
	ClipAtts.plane3Origin = (0, 0, 0)
	ClipAtts.plane1Normal = (1, 0, 0)
	ClipAtts.plane2Normal = (0, 1, 0)
	ClipAtts.plane3Normal = (0, 0, 1)
	ClipAtts.planeInverse = 1
	ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1  # None, Plane1, Plane2, Plane3
	ClipAtts.center = (3.2, 3.55, 3.55)
	ClipAtts.radius = 0.2
	ClipAtts.sphereInverse = 0
	SetOperatorOptions(ClipAtts, 3)

def initializeBloodCell(x,y,z,i):
	OpenDatabase("localhost:/Users/caseychesshir/Dropbox/CIS/CIS410/custom_project/bloodcells/bloodcell copy " + str(i+2) + ".obj", 0)

	AddPlot("Subset", "OBJMesh", 1, 0)
	moveBloodCell(x,y,z)
	#{-7.5, -0.5, -1}, {-5.5, 0.5, 1}
	#{mins}, {maxs}

def moveBloodCell(x,y,z):
	x = float(x)
	y = float(y)
	z = float(z)
	RemoveOperator(1,0)
	AddOperator("Transform", 0)
	TransformAtts = TransformAttributes()
	TransformAtts.doRotate = 0
	TransformAtts.rotateOrigin = (0, 0, 0)
	TransformAtts.rotateAxis = (0, 0, 1)
	TransformAtts.rotateAmount = 0
	TransformAtts.rotateType = TransformAtts.Deg  # Deg, Rad
	TransformAtts.doScale = 0
	TransformAtts.scaleOrigin = (0, 0, 0)
	TransformAtts.scaleX = 10
	TransformAtts.scaleY = 10
	TransformAtts.scaleZ = 10
	TransformAtts.doTranslate = 1
	TransformAtts.translateX = x
	TransformAtts.translateY = y
	TransformAtts.translateZ = z
	TransformAtts.transformType = TransformAtts.Similarity  # Similarity, Coordinate, Linear
	TransformAtts.inputCoordSys = TransformAtts.Cartesian  # Cartesian, Cylindrical, Spherical
	TransformAtts.outputCoordSys = TransformAtts.Spherical  # Cartesian, Cylindrical, Spherical
	TransformAtts.continuousPhi = 0
	TransformAtts.m00 = 1
	TransformAtts.m01 = 0
	TransformAtts.m02 = 0
	TransformAtts.m03 = 0
	TransformAtts.m10 = 0
	TransformAtts.m11 = 1
	TransformAtts.m12 = 0
	TransformAtts.m13 = 0
	TransformAtts.m20 = 0
	TransformAtts.m21 = 0
	TransformAtts.m22 = 1
	TransformAtts.m23 = 0
	TransformAtts.m30 = 0
	TransformAtts.m31 = 0
	TransformAtts.m32 = 0
	TransformAtts.m33 = 1
	TransformAtts.invertLinearTransform = 0
	TransformAtts.vectorTransformMethod = TransformAtts.AsDirection  # None, AsPoint, AsDisplacement, AsDirection
	TransformAtts.transformVectors = 0
	SetOperatorOptions(TransformAtts, 0)
	DrawPlots()

def initializeVectorField():
	AddPlot("Vector","velocity") # why can't this happen in main? 
	VectorAtts = VectorAttributes()
	VectorAtts.glyphLocation = VectorAtts.AdaptsToMeshResolution  # AdaptsToMeshResolution, UniformInSpace
	VectorAtts.useStride = 0
	VectorAtts.stride = 1
	VectorAtts.nVectors = 4000
	VectorAtts.lineStyle = VectorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
	VectorAtts.lineWidth = 0
	VectorAtts.scale = 0.25
	VectorAtts.scaleByMagnitude = 1
	VectorAtts.autoScale = 1
	VectorAtts.headSize = 0.25
	VectorAtts.headOn = 1
	VectorAtts.colorByMag = 1
	VectorAtts.useLegend = 1
	VectorAtts.vectorColor = (0, 0, 0, 255)
	VectorAtts.colorTableName = "Default"
	VectorAtts.invertColorTable = 0
	VectorAtts.vectorOrigin = VectorAtts.Tail  # Head, Middle, Tail
	VectorAtts.minFlag = 0
	VectorAtts.maxFlag = 0
	VectorAtts.limitsMode = VectorAtts.OriginalData  # OriginalData, CurrentPlot
	VectorAtts.min = 0
	VectorAtts.max = 1
	VectorAtts.lineStem = VectorAtts.Line  # Cylinder, Line
	VectorAtts.geometryQuality = VectorAtts.Fast  # Fast, High
	VectorAtts.stemWidth = 0.08
	VectorAtts.origOnly = 1
	VectorAtts.glyphType = VectorAtts.Arrow  # Arrow, Ellipsoid
	SetPlotOptions(VectorAtts)
	DrawPlots()

def processfile(textfile):
	f = open(resultspath + textfile,"r+") # f : file
	initialconditions = []
	datapoints = [[] for i in range(100)]
	numCells = 0
	currentLine = -1 
	for line in f:
		if line[0] == 'S':
			tokens = line.split(" ")
			initialconditions.append([tokens[3],tokens[4],tokens[5]])
			numCells += 1		
			currentLine+=1
		else:
			line.lstrip()
			line.rstrip()
			tokens = line.split(" ")
			datapoints[currentLine].append(float(tokens[0]))
			datapoints[currentLine].append(float(tokens[1]))
			datapoints[currentLine].append(float(tokens[2]))
	
	counter = 0
	for point in initialconditions[:30]:
		initializeBloodCell(point[0],point[1],point[2],counter)		
		counter+=1
	
	'''
	for j in range(0,len(datapoints[i]),3):
		for i in range(numCells):
			moveBloodCell(datapoints[j][i], datapoints[j][i+1], datapoints[j][i+2])
	
	for i in datapoints:
		fly(i)
	'''
	fly(datapoints[0])
	for step in range(0,100,3):
		for i in range(numCells):
			pass
			#initializeBloodCell(datapoints[i][step],datapoints[i][step+1],datapoints[i][step+2],i)

def cross(a,b):
	c = [a[1]*b[2] - a[2]*b[1],
		 a[2]*b[0] - a[0]*b[2],
		 a[0]*b[1] - a[1]*b[0]]
	return c

def fly(points):
	a = points[0]
	b = points[1]
	g = points[2]

	d = points[3]
	e = points[4]
	f = points[5]


	j = d-a
	k = e-b
	l = f-g

	x = []
	cpts = [0 for i in range(len(points)/3)]
	for i in range(len(points)/3):
		x = x + [float(i) / float(len(points))]
		cpts[i] = View3DAttributes()
		#cpts[i].viewNormal = (points[3*i],points[3*i+1],points[3*i+2])
		#cpts[i].focus = (4,4,5)
		#cpts[i].viewNormal = (4,4,5)
		#cpts[i].viewNormal = (points[3*i-3],points[3*i-2],points[3*i-1])
		cpts[i].viewNormal = (0.6790127871442293, -0.7229199320757748, -0.12778265415220738)

		cpts[i].focus = (points[3*i],points[3*i+1],points[3*i+2])
		#cpts[i].viewUp = (c1,c2,c3)
		cpts[i].viewUp = (1,1,1)
		cpts[i].viewAngle=30
		cpts[i].parallelScale = 1
		cpts[i].nearPlane = -4
		cpts[i].farPlane = 4
		cpts[i].perspective = 1
	#j = j/(d-a)
	#k = k/(e-b)
	#l = l/(f-g)
	cpts[-1] = cpts[0]

	nsteps = 1000
	for i in range(nsteps/2-1):
		t = float(i) / float(nsteps - 1)
		len(cpts)
		len(x)
		try:
			c = EvalCubicSpline(t, x, cpts)
			c.imageZoom = 3
			moveBloodCell(c.focus[0],c.focus[1],c.focus[2])
			SetView3D(c)
		except Exception as e:
			print(e)

		

def main():
	Source("Dropbox/CIS/CIS410/custom_project/initialize.py") #have to resource this file for some reason 
	DeleteAllPlots()
	ResetView()
	initializeMesh()
	#initializeStreamline()
	#initializeVectorField() # this can't happen in main for some reason 
	initializeClip1()
	initializeClip2()
	initializeClip3()
	start = positions[0:3] # first 3D location 

	#for i in range(0,1000,3):
	#	moveBloodCell(positions[i],positions[i+1],positions[i+2])
	#	DrawPlots()
	#	print(i)

	OpenDatabase("localhost:/Users/caseychesshir/Dropbox/CIS/CIS410/custom_project/bloodcell.obj", 0) # have to reload for some reason 
	processfile("bigresults.txt")	
	
	initializeStreamline()
	DrawPlots()


