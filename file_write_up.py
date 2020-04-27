
import hashlib
import os
import time
import random


#This class(??) does the writing of the vxa file that will be read by voxelyze and used in the testing!
#  I think at some point I will have to work out how to put in indiviual bitsof the phenotype in, but we'll leave it for the moment.  

def write_voxelyze_file(t,morphology,stiffness):
  
    bump1=0.1+(0.1*t)
    bump2=15+(5*t)

    #this first part is making sure that our morphology and stiffness data is in the correct 
    #format to be read by voxelyze. It also makes sure that the muscle core remains constant and its sculpted away
    t=str(t)
    morphology_l1=str(morphology[0])
    morphology_l1=morphology_l1.replace(',','')
    morphology_l1=morphology_l1.replace(" ", "")

    morphology_l2=morphology[1]

    for i in range(12,24):
        morphology_l2[i]=4


    morphology_l2=str(morphology_l2)
    morphology_l2 = morphology_l2.replace(',','')
    morphology_l2=morphology_l2.replace(" ", "")
    
    
    morphology_l3=morphology[2]
    for i in range(12,24):
        morphology_l3[i]=4


    morphology_l3=str(morphology_l3)
    morphology_l3 = morphology_l3.replace(',','')
    morphology_l3=morphology_l3.replace(" ", "")

    morphology_l4=morphology[3]
    for i in range(12,24):
        morphology_l4[i]=4


            
    morphology_l4=str(morphology_l4)
    morphology_l4 = morphology_l4.replace(',','')
    morphology_l4=morphology_l4.replace(" ", "")
    

    morphology_l5=morphology[4]
    for i in range(12,24):
        morphology_l5[i]=4



    morphology_l5=str(morphology_l5)
    morphology_l5 = morphology_l5.replace(',','')
    morphology_l5=morphology_l5.replace(" ", "")
    

    morphology_l6=str(morphology[5])
    morphology_l6 = morphology_l6.replace(',','')
    morphology_l6=morphology_l6.replace(" ", "")


    stiffness_l1=str(stiffness[0])
    stiffness_l2=stiffness[1]
    for i in range(12,24):
        stiffness_l2[i]=500000

    stiffness_l3=stiffness[2]
    for i in range(12,24):
        stiffness_l3[i]=500000

    stiffness_l2=str(stiffness_l2)
    stiffness_l4=stiffness[3]
    for i in range(12,24):
        stiffness_l4[i]=500000

    
    stiffness_l3=str(stiffness_l3)
    stiffness_l4=str(stiffness_l4)
    stiffness_l5=stiffness[4]
    for i in range(12,24):
        stiffness_l5[i]=500000

    stiffness_l5=str(stiffness_l5)	
    stiffness_l6=str(stiffness[5])

    #first we need to open the file, for writing purposes
    voxelyze_file = open("katt"+t+".vxa", "w")

    #then this is the actual writing bit
    voxelyze_file.write(
    "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n\
    <VXA Version=\"1.0\">\n\
    <Simulator>\n\
    <Integration>\n\
    <Integrator>0</Integrator>\n\
    <DtFrac>0.9</DtFrac>\n\
    </Integration>\n\
    <Damping>\n\
    <BondDampingZ>1</BondDampingZ>\n\
    <ColDampingZ>0.8</ColDampingZ>\n\
    <SlowDampingZ>0.01</SlowDampingZ>\n\
    </Damping>\n\
    <Collisions>\n\
    <SelfColEnabled>1</SelfColEnabled>\n\
    <ColSystem>3</ColSystem>\n\
    <CollisionHorizon>2</CollisionHorizon>\n\
    </Collisions>\n\
    <Features>\n\
    <FluidDampEnabled>0</FluidDampEnabled>\n\
    <PoissonKickBackEnabled>0</PoissonKickBackEnabled>\n\
    <EnforceLatticeEnabled>0</EnforceLatticeEnabled>\n\
    </Features>\n\
    <SurfMesh>\n\
    <CMesh>\n\
    <DrawSmooth>1</DrawSmooth>\n\
    <Vertices/>\n\
    <Facets/>\n\
    <Lines/>\n\
    </CMesh>\n\
    </SurfMesh>\n\
    <StopCondition>\n\
    <StopConditionType>2</StopConditionType>\n\
    <StopConditionValue>8</StopConditionValue>\n\
    <InitCmTime>1.0</InitCmTime>\n\
    </StopCondition>\n\
    <GA>\n\
    <WriteFitnessFile>1</WriteFitnessFile>\n\
    <FitnessFileName>my_fitness"+t+".xml</FitnessFileName>\n\
    <QhullTmpFile>Qhull_temp"+t+"</QhullTmpFile>\n\
    <CurvaturesTmpFile>curve_temp"+t+"</CurvaturesTmpFile>\n\
    </GA>\n\
    </Simulator>\n\
    <Environment>\n\
    <Fixed_Regions>\n\
    <NumFixed>0</NumFixed>\n\
    </Fixed_Regions>\n\
    <Forced_Regions>\n\
    <NumForced>0</NumForced>\n\
    </Forced_Regions>\n\
    <Gravity>\n\
    <GravEnabled>1</GravEnabled>\n\
    <GravAcc>-9.81</GravAcc>\n\
    <FloorEnabled>1</FloorEnabled>\n\
    <FloorSlope>0</FloorSlope>\n\
    <bump_size>bump1</bump_size>\n\
    <bump_seperation>bump2</bump_seperation>\n\
    </Gravity>\n\
    <Thermal>\n\
    <TempEnabled>1</TempEnabled>\n\
    <TempAmp>39</TempAmp>\n\
    <TempBase>25</TempBase>\n\
    <VaryTempEnabled>1</VaryTempEnabled>\n\
    <TempPeriod>0.25</TempPeriod>\n\
    </Thermal>\n\
    </Environment>\n\
    <VXC Version=\"0.93\">\n\
    <Lattice>\n\
    <Lattice_Dim>0.05</Lattice_Dim>\n\
    <X_Dim_Adj>1</X_Dim_Adj>\n\
    <Y_Dim_Adj>1</Y_Dim_Adj>\n\
    <Z_Dim_Adj>1</Z_Dim_Adj>\n\
    <X_Line_Offset>0</X_Line_Offset>\n\
    <Y_Line_Offset>0</Y_Line_Offset>\n\
    <X_Layer_Offset>0</X_Layer_Offset>\n\
    <Y_Layer_Offset>0</Y_Layer_Offset>\n\
    </Lattice>\n\
    <Voxel>\n\
    <Vox_Name>BOX</Vox_Name>\n\
    <X_Squeeze>1</X_Squeeze>\n\
    <Y_Squeeze>1</Y_Squeeze>\n\
    <Z_Squeeze>1</Z_Squeeze>\n\
    </Voxel>\n\
    <Palette>\n\
    <Material ID=\"1\">\n\
    <MatType>0</MatType>\n\
    <Name>Passive_Soft</Name>\n\
    <Display>\n\
    <Red>0</Red>\n\
    <Green>1</Green>\n\
    <Blue>1</Blue>\n\
    <Alpha>1</Alpha>\n\
    </Display>\n\
    <Mechanical>\n\
    <MatModel>0</MatModel>\n\
    <Elastic_Mod>1000</Elastic_Mod>\n\
    <Plastic_Mod>0</Plastic_Mod>\n\
    <Yield_Stress>0</Yield_Stress>\n\
    <FailModel>0</FailModel>\n\
    <Fail_Stress>0</Fail_Stress>\n\
    <Fail_Strain>0</Fail_Strain>\n\
    <Density>1200.0</Density>\n\
    <Poissons_Ratio>0.4</Poissons_Ratio>\n\
    <CTE>0</CTE>\n\
    <uStatic>1</uStatic>\n\
    <uDynamic>0.5</uDynamic>\n\
    </Mechanical>\n\
    </Material>\n\
    <Material ID=\"2\">\n\
    <MatType>0</MatType>\n\
    <Name>Passive_Hard</Name>\n\
    <Display>\n\
    <Red>0</Red>\n\
    <Green>0</Green>\n\
    <Blue>1</Blue>\n\
    <Alpha>1</Alpha>\n\
    </Display>\n\
    <Mechanical>\n\
    <MatModel>0</MatModel>\n\
    <Elastic_Mod>10000000</Elastic_Mod>\n\
    <Plastic_Mod>0</Plastic_Mod>\n\
    <Yield_Stress>0</Yield_Stress>\n\
    <FailModel>0</FailModel>\n\
    <Fail_Stress>0</Fail_Stress>\n\
    <Fail_Strain>0</Fail_Strain>\n\
    <Density>2200.0</Density>\n\
    <Poissons_Ratio>0.4</Poissons_Ratio>\n\
    <CTE>0</CTE>\n\
    <uStatic>1</uStatic>\n\
    <uDynamic>0.5</uDynamic>\n\
    </Mechanical>\n\
    </Material>\n\
    <Material ID=\"3\">\n\
    <MatType>0</MatType>\n\
    <Name>Active_+</Name>\n\
    <Display>\n\
    <Red>1</Red>\n\
    <Green>0</Green>\n\
    <Blue>0</Blue>\n\
    <Alpha>1</Alpha>\n\
    </Display>\n\
    <Mechanical>\n\
    <MatModel>0</MatModel>\n\
    <Elastic_Mod>1.0e+006</Elastic_Mod>\n\
    <Plastic_Mod>0</Plastic_Mod>\n\
    <Yield_Stress>0</Yield_Stress>\n\
    <FailModel>0</FailModel>\n\
    <Fail_Stress>10</Fail_Stress>\n\
    <Fail_Strain>0</Fail_Strain>\n\
    <Density>1200.0</Density>\n\
    <Poissons_Ratio>0.4</Poissons_Ratio>\n\
    <CTE>0</CTE>\n\
    <uStatic>1</uStatic>\n\
    <uDynamic>0.5</uDynamic>\n\
    </Mechanical>\n\
    </Material>\n\
    <Material ID=\"4\">\n\
    <MatType>0</MatType>\n\
    <Name>Active_-</Name>\n\
    <Display>\n\
    <Red>0</Red>\n\
    <Green>1</Green>\n\
    <Blue>0</Blue>\n\
    <Alpha>1</Alpha>\n\
    </Display>\n\
    <Mechanical>\n\
    <MatModel>0</MatModel>\n\
    <Elastic_Mod>1.0e+006</Elastic_Mod>\n\
    <Plastic_Mod>0</Plastic_Mod>\n\
    <Yield_Stress>0</Yield_Stress>\n\
    <FailModel>0</FailModel>\n\
    <Fail_Stress>0</Fail_Stress>\n\
    <Fail_Strain>0</Fail_Strain>\n\
    <Density>1200.0</Density>\n\
    <Poissons_Ratio>0.4</Poissons_Ratio>\n\
    <CTE>0.02</CTE>\n\
    <uStatic>1</uStatic>\n\
    <uDynamic>0.5</uDynamic>\n\
    </Mechanical>\n\
    </Material>\n\
    <Material ID=\"5\">\n\
    <MatType>0</MatType>\n\
    <Name>Aperture</Name>\n\
    <Display>\n\
    <Red>1</Red>\n\
    <Green>0.784</Green>\n\
    <Blue>0</Blue>\n\
    <Alpha>1</Alpha>\n\
    </Display>\n\
    <Mechanical>\n\
    <MatModel>0</MatModel>\n\
    <Elastic_Mod>5e+007</Elastic_Mod>\n\
    <Plastic_Mod>0</Plastic_Mod>\n\
    <Yield_Stress>0</Yield_Stress>\n\
    <FailModel>0</FailModel>\n\
    <Fail_Stress>0</Fail_Stress>\n\
    <Fail_Strain>0</Fail_Strain>\n\
    <Density>1200.0</Density>\n\
    <Poissons_Ratio>0.4</Poissons_Ratio>\n\
    <CTE>-0.04</CTE>\n\
    <uStatic>1</uStatic>\n\
    <uDynamic>0.5</uDynamic>\n\
    </Mechanical>\n\
    </Material>\n\
    </Palette>\n\
    <Structure Compression=\"ASCII_READABLE\">\n\
    <X_Voxels>6</X_Voxels>\n\
    <Y_Voxels>6</Y_Voxels>\n\
    <Z_Voxels>6</Z_Voxels>\n\
    <Data>\n\
    <Layer><![CDATA"+morphology_l1+"]></Layer>\n\
    <Layer><![CDATA"+morphology_l2+"]></Layer>\n\
    <Layer><![CDATA"+morphology_l3+"]></Layer>\n\
    <Layer><![CDATA"+morphology_l4+"]></Layer>\n\
    <Layer><![CDATA"+morphology_l5+"]></Layer>\n\
    <Layer><![CDATA"+morphology_l6+"]></Layer>\n\
    </Data>\n\
    <PhaseOffset>\n\
    <Layer><![CDATA[0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,]]></Layer>\n\
    <Layer><![CDATA[0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,]]></Layer>\n\
    <Layer><![CDATA[0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,]]></Layer>\n\
    <Layer><![CDATA[0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,]]></Layer>\n\
    <Layer><![CDATA[0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,]]></Layer>\n\
    <Layer><![CDATA[0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,0,-0.1,-0.2,-0.3,-0.4,-0.5,]]></Layer>\n\
    </PhaseOffset>\n\
    <Stiffness>\n\
    <MinElasticMod>1000.0</MinElasticMod>\n\
    <MaxElasticMod>1000000</MaxElasticMod>\n\
    <Layer><![CDATA"+stiffness_l1+"]></Layer>\n\
    <Layer><![CDATA"+stiffness_l2+"]></Layer>\n\
    <Layer><![CDATA"+stiffness_l3+"]></Layer>\n\
    <Layer><![CDATA"+stiffness_l4+"]></Layer>\n\
    <Layer><![CDATA"+stiffness_l5+"]></Layer>\n\
    <Layer><![CDATA"+stiffness_l6+"]></Layer>\n\
    </Stiffness>\n\
    </Structure>\n\
    </VXC>\n\
    </VXA>")
    voxelyze_file.close()
    return

def read_voxelyze_file(ta):
    this_file=open("my_fitness"+ta+".xml")
    tag="<normDistY>"

    results=0
    results_X=0
    results_x=0
    for line in this_file:
                if tag in line:
                    results_Y = float(line[line.find(tag) + len(tag):line.find("</" + tag[1:])])
                    results=abs(results_Y)
    this_file.close()

    this_file.close()
    return results
def read_voxelyze_file_un(ta):
    this_file=open("my_fitness"+ta+".xml")
    tag="<normDistY>"

    results_Y=0
    results_X=0
    results_x=0
    for line in this_file:
                if tag in line:
                    results_Y = float(line[line.find(tag) + len(tag):line.find("</" + tag[1:])])
                    results_Y=abs(results_Y)
    this_file.close()
    this_file=open("my_fitness"+ta+".xml")
    tag="<normDistX>"
    for line in this_file:
                if tag in line:
                    results_X = float(line[line.find(tag) + len(tag):line.find("</" + tag[1:])])
                    results_X=(results_X)            

    results=[results_Y,results_X]

    this_file.close()
   
    return results  
  
