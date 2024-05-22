#pragma rtGlobals=3		// Use modern global access method and strict wave access.
///////////////////////////////////////////////////////////////////////////////////
//////////// Guilherme Neves version 06/24
////// e-mail: guilherme.neves@kcl.ac.uk
////	Juan Burrone Lab
//// MRC Center for NeuroDevelopmental Disorders
/// IoPPN, King's College London
//// Compatible with IgorPro versions 6-9

/////////////////////////////////////////////////////////////////////////////////////
//// Menus to call functions from Macros Menu
Menu "Macros"
	SubMenu "Load Data"
		"Analyse all files on folder", AnalyseAllBranches()
		"Normalize Fluor per cell", NormalizeFluoMedFolders()
		"Check Results of Branch Analysis", Diagnostics()
	End
	SubMenu "Analysis"
		"Summarize Branch Parameters", Summarize_Sel_Fold()
	End
		
End
//////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////Folder making Igor///////////////////////////

/////////// Load all files to experiment
// Called from the Macros Menu using the "Analyse all files on folder" tab
///////// Under Load Data SubMenu
//// Example File structure in Example Dataset Folder
//// Data Obtained using Fiji functions SNT (for swc files)
/// and ROI Manager MultiMeasure (for csv files)	
/// Calls function Analyze_Puncta_Folder
/////////////////////////////////////////////////////////////////////////
Function AnalyseAllBranches()

Variable Zscale=0.5
Prompt Zscale, "Z step of images - microns"
DoPrompt "Input Z scale", Zscale
		
SetDataFolder root:
String pathName = "TestPath"
NewPath/Q/M="Choose Folder containing Analysis files" $pathName

Analyze_Puncta_Folder(pathName,Zscale)

End
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// Load all files to experiment
//// Example File structure in Example Dataset Folder
//// Data Obtained using Fiji functions SNT (for swc files)
/// and ROI Manager MultiMeasure (for csv files)	
/// Calls functions LoadallSwcFromFolder (for Loading swc files)
/// and LoadWaves_Synapses (for Loading csv files)
///	and CreateIgorFolders (for Creating Folders in Igor's Data Browser)
Function Analyze_Puncta_Folder(path,Z_scale)
String path
Variable Z_scale

//// Swc File Loading (all files ending in swc in folder loaded)
LoadallSwcFromFolder(path)

String list = IndexedFile($path, -1, "????"), nextFilename, branch_name, nameExc, nameInh
Variable i=0

//////// Swc File Loading (all files ending in csv in folder loaded)
////// Files ending in tdt.csv are inhibitory synapse measurements
////// Files ending in gfp.csv are excitatory synapse measurements
list = SortList(list, ";", 16)
do
	nextFilename=StringFromList(i, list)
	if(strlen(nextFilename)==0)
		break
	endif
	if (stringmatch(nextFilename,"*.csv")==1)
		if (stringmatch(nextFilename,"*tdt*")==1)
			branch_name=removeending(nextFilename,"tdt.csv")
			nameInh=branch_name+"_Inh_"
			Wave wTree=$(branch_name+"_swc")		
			LoadWaves_Synapses(nameInh,Z_scale,wTree,path,nextFilename)
		elseif (stringmatch(nextFilename,"*gfp*")==1)
			branch_name=removeending(nextFilename,"gfp.csv")
			nameExc=branch_name+"_Exc_"
			Wave wTree=$(branch_name+"_swc")	
			LoadWaves_Synapses(nameExc,Z_scale,wTree,path,nextFilename)
		endif
	endif
	i+=1		
	while(1)
/// Creating Folders in Igor's Data Browser (one Folder per Branch)
/// Note that Cell File name structures of swc and csv files need to match
CreateIgorFolders()

End
//////////////////////////////////////		Load SWC FIles	///////////////////////////////////////////////////////////////
//////////// /// Calls function Load_SWC_SNT_Folder to load each swc file in folder
Function LoadallSwcFromFolder(path)
String path

String list = IndexedFile($path, -1, ".swc"), nextFilename
Variable i=0

list = SortList(list, ";", 16)
do
	nextFilename=StringFromList(i, list)
	if(strlen(nextFilename)==0)
		break
	endif
	Load_SWC_SNT_Folder(path,nextFilename)
	i+=1		
while(1)
	
End
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////// For each swc file, loads 3D co-ordinates and measure euclidian distance
//////// to branch origin following the branch path using function DistanceBtweenNodes
Function Load_SWC_SNT_Folder(path,file)
String path, file

String name="test"
LoadWave/J/M/D/N=$name/O/Q/V={"\t, "," $",0,0}/K=1/L={0,6,0,2,3}/P=$path file

String namefile=Removeending(file,"-000.swc")
Rename test0, $(namefile+"_swc")
Wave wSWC=$(namefile+"_swc")
String labelX="X_coord", labelY="Y_coord", labelZ="Z_coord"
SetDimLabel 1, 0, $labelX, wSWC
SetDimLabel 1, 1, $labelY, wSWC
SetDimLabel 1, 2, $labelZ, wSWC
DistanceBtweenNodes(wSWC)
End

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// For each swc file, this function measures 3D euclidian distance
//////// to branch origin following the branch path for each branch node
/// Uses Function DistanceXYZ
Function DistanceBtweenNodes(w1)
Wave w1

Variable nColumns=DimSize(w1,1), nrows=DimSize(w1,0),i
Variable x1,x2,y1,y2,z1,z2,D=0
String labelDist="Dist_Origin"
InsertPoints/M=1 nColumns, 1,w1
w1[0][nColumns]=0
SetDimLabel 1, nColumns, $labelDist, w1
For (i=1;i<nrows;i+=1)
	x1=w1[i-1][0]
	y1=w1[i-1][1]
	z1=w1[i-1][2]
	x2=w1[i][0]
	y2=w1[i][1]
	z2=w1[i][2]
	D+=DistanceXYZ(x1,y1,z1,x2,y2,z2)
	w1[i][nColumns]=D
EndFor

End
////////////////////////////////////////////	Ends SWC Loading 		//////////////////////////////////////////////////////////////////
///////////////////		 Loads Synapse Measurement File			//////////////////////////////////////////////////////////////
//////////// 
//////////// Loads csv files. Matches their 3D location with the closest branch
//////////// Node (using minimum 3D euclidian distance), and assigns a distance of the synapse to the branch origin
/////////// Uses Functions DistanceXY and DistanceXYZ
////////// Also identifies Z position of ROI by finding the maximum fluorescence in adjacent
////////// sections
Function LoadWaves_Synapses(name,Z_scale,wTree,path,File)
Wave wTree
Variable Z_scale
String name, File,path

//Load Results Table
//Set Measurements in ImageJ are for Min and Max gray value and Centroid
//Measurements using Filtered Median3D

//For Fiji Files that generate Excell formatted results use
LoadWave/Q/A/J/D/K=1/L={0,1,0,1,0}/P=$path File
Variable nwaves=V_flag, nROIS=(nwaves)/4
String ListofWaves=S_waveNames

// Make Table to extract values into columns (0-XROI; 1-YROI; 2-ZROI; 3-index; 4-MaxFluo; 5-SliceROI)
// (6-TreeDistance; 7 - MinDistancetoROI; 8-NodeX; 9-NodeY; 10-NodeZ; 11-NodeSl; 12-NodeAppr
Make/O/N=(nROIS,13) $(name+"Results")
Wave wResults=$(name+"Results")
String l_XROI="X_centroid", l_YROI="Y_centroid",l_ZROI="Z_centroid", l_iROI="index_ROI", l_MaxFluo="Max_Fluo"
String l_SliceROI="Slice_ROI", l_Dist="Distance to Soma",l_MinDist="Dist_ROI_node"
String l_Xnd="X_node", l_Ynd="Y_node",l_Znd="Z_node", l_Slicend="Slice_node", l_SlicendAppr="Z_nodeXYmatch"
SetDimLabel 1, 0, $l_XROI, wResults
SetDimLabel 1, 1, $l_YROI, wResults
SetDimLabel 1, 2, $l_ZROI, wResults
SetDimLabel 1, 3, $l_iROI, wResults
SetDimLabel 1, 4, $l_MaxFluo, wResults
SetDimLabel 1, 5, $l_SliceROI, wResults
SetDimLabel 1, 6, $l_Dist, wResults
SetDimLabel 1, 7, $l_MinDist, wResults
SetDimLabel 1, 8, $l_Xnd, wResults
SetDimLabel 1, 9, $l_Ynd, wResults
SetDimLabel 1, 10, $l_Znd, wResults
SetDimLabel 1, 11, $l_Slicend, wResults
SetDimLabel 1, 12, $l_SlicendAppr, wResults

// Extract X and Y co-ordinates from Excell File
String theWave
Variable index=0, index_wave=0, ROI_index=0

For (ROI_index=0; ROI_index<nROIS; ROI_index+=1)
	wResults[ROI_index][3]=ROI_index+1
	For (index_wave=0; index_wave<4; index_wave+=1)
		index=ROI_index*4+index_wave	
		switch(index_wave)			
			case 0:
			// Get the next wave name
			//Delete wave with Minimum
			theWave = StringFromList(index, ListofWaves)
			Wave waveminim=$theWave
			KillWaves waveminim
			break
			// Ignore Maximum wave for now
			case 2:				
			// Get the next wave name
			// Get X co-ordinate
			theWave = StringFromList(index, ListofWaves)
			Wave waveX=$theWave
			wResults[ROI_index][0]=waveX[0]
			KillWaves waveX
			break
			case 3:				
			// Get the next wave name
			// Get Y co-ordinate
			theWave = StringFromList(index, ListofWaves)
			Wave waveY=$theWave
			wResults[ROI_index][1]=waveY[0]
			KillWaves waveY
			break
		endswitch	
	EndFor	
EndFor

// Get the Approx Z for ROI from node position (use XY distance only)
Variable number_nodes=DimSize(wTree,0)
Variable Znode, index1, x1, y1, x2, y2, Zappr
Make/O/N=(number_nodes) $"Dist"
Wave wdist=$"Dist"

// Find Minimal Distance ROI node in XY co-ordinates
For (index=0; index<NROIS; index+=1)
	x2=wResults[index][0]
	y2=wResults[index][1]
	For (index1=0; index1<number_nodes; index1+=1)
		x1=wTree[index1][0]
		y1=wTree[index1][1]
		wdist[index1]=DistanceXY(x1,y1,x2,y2)
	EndFor
	WaveStats/Q wdist
	Zappr=wTree[V_minloc][2]
	Znode=round(Zappr/Z_scale)
	wResults[index][12]=Znode
EndFor
/////////////////////////////////////////////		Important Variable	//////////////////////////////////////////////////////////////
// Get the maximum Fluorescence and Z position from the 12 slices (4 microns) adjacent to the node
// Can be modified using variable VZsearch = number of slices to search for maximum around the node

	Variable VZsearch=12
/////////////////////////////////////////////////////////////////////////////////////////////////////
// Get Max fluorescence wave and extract Max and Z position of ROI
For (ROI_index=0; ROI_index<nROIS; ROI_index+=1)
	Znode=wResults[ROI_index][12]
	For (index_wave=0; index_wave<4; index_wave+=1)
		index=ROI_index*4+index_wave	
		switch(index_wave)			
			case 1:
			// Get the next wave name
			theWave = StringFromList(index, ListofWaves)
			Wave wavemaxim=$theWave
			WaveStats/Q/R=[(Znode-VZsearch),(Znode+VZsearch) ] wavemaxim
			wResults[ROI_index][4]=V_max
			wResults[ROI_index][5]=V_maxloc
			wResults[ROI_index][2]=(V_maxloc)*Z_scale
			KillWaves wavemaxim
			break
		endswitch	
	EndFor	
EndFor

// Correlates ROI coordinates with Tree Nodes in 3D
Variable z1,z2
For (index=0; index<nROIS; index+=1)
	// This Part creates a wave with the distance from the ROI to all the nodes
	x2=wResults[index][0]
	y2=wResults[index][1]
	z2=wResults[index][2]
	 // runs through nodes
	For (index1=0; index1<number_nodes; index1+=1)
		x1=wTree[index1][0]
		y1=wTree[index1][1]
		z1=wTree[index1][2]
		wdist[index1]=DistanceXYZ(x1,y1,z1,x2,y2,z2)		
	EndFor
	// Wave wdist created
	WaveStats/Q wdist
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//To prevent crashes when coordinates are incorrect //before next four lines
	if (V_min>5)
		String ErrorMess="On "+File+". Wrong co-ordinates ROI too far from the tree"
		Abort ErrorMess
	endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	wResults[index][7]=V_min
	// Minimum distance to node recorded for quality checking
	// Wave name MinD
	// Location of selected node (X, Y and Z in microns recorded)
	// Slice containing the closest node (Wave ending in Sl)
	// Distance from node to beggining of trace recorded (Wave ending with Dist)
	wResults[index][6]=wTree[V_minloc][3]
	wResults[index][8]=wTree[V_minloc][0]
	wResults[index][9]=wTree[V_minloc][1]
	wResults[index][10]=wTree[V_minloc][2]
	wResults[index][11]=round(wResults[index][10]/Z_scale)
EndFor
KillWaves wdist

End
/////////////////////////////////		End Load Measurement Results	//////////////////////////////////////////////////
//////////////		Organize Data Into Folders	in Igor Data Browser		////////////////////////////////////////////////
Function CreateIgorFolders()

String RootfldrSav= GetDataFolder(1)
String list_swc=Wavelist("*_swc",";",""), list_Exc=Wavelist("*Exc*",";",""), list_Inh=Wavelist("*Inh*",";","")
String Branchname, swcname, Inhname, Excname
Variable i=0
list_swc = SortList(list_swc, ";", 16)
list_Exc =  SortList(list_Exc, ";", 16)
list_Inh =  SortList(list_Inh, ";", 16)
	do
		swcname=StringFromList(i, list_swc)
		if(strlen(swcname)==0)
			break
		endif
		Branchname=removeending(swcname,"_swc")
		Excname=StringFromList(i, list_Exc)
		Inhname=StringFromList(i, list_Inh)
		Wave wswc=$swcname
		Wave wExc=$Excname
		Wave wInh=$Inhname
		NewDataFolder/O $Branchname
		DFREF dfr = root:$Branchname
		MoveWave wswc, dfr
		MoveWave wExc, dfr
		MoveWave wInh, dfr
		i+=1		
	while(1)

End
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////		Find 2D Euclidian distance between 2 points			////////////////////////////////////////////////

Function DistanceXY(x1,y1,x2,y2)
Variable x1,y1,x2,y2

Variable dist=sqrt((x1-x2)^2+(y1-y2)^2)
return dist
End
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////		Find 3D Euclidian distance between 2 points			////////////////////////////////////////////////

Function DistanceXYZ(x1,y1,z1,x2,y2,z2)
Variable x1,y1,z1,x2,y2,z2

Variable dist=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
return dist
End
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////		Normalize Fluorescence across image			////////////////////////////////////////////////
////////// Called from the Macros Menu using the "Normalize Fluor per cell" tab
///////// Under Load Data SubMenu
/////////// Normalizes the fluorescence intensity of all ROIs in image of the same  /////////
///// type (Excitation or inhibition) to the median of the distribution
////////// Uses Function GetBrowserSelection_S and AutExtractWaveFromTable
/////////// CreateResultsTableInh and CreateResultsTableExc
////////// All folders containing measurements from the same cell
///////// Need to be selected in Igor's Data Browser
//////// Creates Summary Tables for Excitation and Inhibition
///////// To be loaded into MATLAB
//////// Summary Tables contain one row per measured synapse
/////// Columns are Distance to branch origin (in microns)
/////// and Normalized Fluorescence
///////// Inhibitory Tables also contain a third column with total branch length
Function NormalizeFluoMedFolders()

SetDataFolder root:
String FolderList=GetBrowserSelection_S(""), dataFolderPathStr, nameFolder, ExcFlNm="FluoExcNm", InhFlNm="FluoInhNm"
//print FolderList
Make/O/N=1 $"NormFlExc", $"NormFlInh"
Wave wNExc=$"NormFlExc", wNInh=$"NormFlInh"
Variable i=0, MaxFluo
String regexp=":([[:alnum:]]+):$"
do
		dataFolderPathStr = StringFromList(i, FolderList, ";")
		if (strlen(dataFolderPathStr) == 0)
			break													// No more data folders.
		endif
		SplitString/E= regexp dataFolderPathStr, nameFolder
		SetDataFolder dataFolderPathStr
		Wave wExc=$(nameFolder+"_Exc_Results")
		Wave wInh=$(nameFolder+"_Inh_Results")
		
		Wave wExcFluoNorm=AutExtractWaveFromTable(wExc,4)
		Rename tempwave, $ExcFlNm 
		MaxFluo=WaveMax(wExcFluoNorm)
//		if (MaxFluo<256)
//			Abort "Excitatory Values extracted from 8 bit image"
//		endif
		Wave wInhFluoNorm=AutExtractWaveFromTable(wInh,4)
		Rename tempwave, $InhFlNm
		MaxFluo=WaveMax(wInhFluoNorm)
//		if (MaxFluo<256)
//			Abort "Inhibitory Values extracted from 8 bit image"
//		endif
		MoveWave wExcFluoNorm, root:
		MoveWave wInhFluoNorm, root:
		SetDataFolder root:
		Concatenate/NP/KILL {wExcFluoNorm}, wNExc
		Concatenate/NP/KILL {wInhFluoNorm}, wNInh
		i += 1		
while(1)
DeletePoints 0,1, wNExc
DeletePoints 0,1, wNInh
Variable MedFluoExc=StatsMedian(wNExc)
Variable MedFluoInh=StatsMedian(wNInh)
KillWaves wNExc, wNInh
i=0
do
		dataFolderPathStr = StringFromList(i, FolderList, ";")
		if (strlen(dataFolderPathStr) == 0)
			break													// No more data folders.
		endif
	
		SplitString/E=regexp dataFolderPathStr, nameFolder
		SetDataFolder $dataFolderPathStr
		Wave wExc=$(nameFolder+"_Exc_Results")		
		Wave wExcFluoNorm=AutExtractWaveFromTable(wExc,4)
		Rename tempwave, $ExcFlNm
		wExcFluoNorm/=MedFluoExc
		if (DimSize(wExc,1)==13)
			InsertPoints/M=1 13,1, wExc
			wExc[][13]=wExcFluoNorm[p]
		elseif (DimSize(wExc,1)==14)
			wExc[][13]=wExcFluoNorm[p]
		else	
			Abort "Wrong Number of Columns in Exc Table"
		endif
		Wave wInh=$(nameFolder+"_Inh_Results")				
		Wave wInhFluoNorm=AutExtractWaveFromTable(wInh,4)
		Rename tempwave, $InhFlNm
		wInhFluoNorm/=MedFluoInh
		if (DimSize(wInh,1)==13)
			InsertPoints/M=1 13,1, wInh
			wInh[][13]=wInhFluoNorm[p]
		elseif (DimSize(wInh,1)==14)
			wInh[][13]=wInhFluoNorm[p]
		else
			Abort "Wrong Number of Columns in Inh Table"
		endif
		Wave wTree=$(nameFolder+"_swc")
		CreateResultsTableInh(wTree,wInh)
		CreateResultsTableExc(wExc)
		KillWaves wInhFluoNorm, wExcFluoNorm
		i+=1
while(1)

SetDataFolder root:

End
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function/WAVE CreateResultsTableExc(wResults)
Wave wResults

Variable nROIs=DimSize(wResults,0)

Make/O/N=(NROIs,2) $("Tbl_Exc")
Wave wSum=$("Tbl_Exc")
wSum[][0]=wResults[p][6]
wSum[][1]=wResults[p][13]
return wSum
End
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function/WAVE CreateResultsTableInh(wTree,wResults)
Wave wTree,wResults

Wave wdist=AutExtractWaveFromTable(wTree,3)
Variable nROIs=DimSize(wResults,0), length=WaveMax(wdist)

Make/O/N=(NROIs,3) $("Tbl_Inh")
Wave wSum=$("Tbl_Inh")
wSum[0][2]=length
wSum[][0]=wResults[p][6]
wSum[][1]=wResults[p][13]
KillWaves wdist
return wSum
End


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GetBrowserSelection_S Function
// Gets the Folder selections from the data Browser 
Function/S GetBrowserSelection_S(WavePathList)
String WavePathList
	
If (strlen(WavePathList) == 0)
	Variable i=0
	String nextPath
	do
		nextPath=GetBrowserSelection(i)
		if(strlen(nextPath)==0)
			break
		endif
		WavePathList+=nextPath+";"
		i+=1
	while(1)
EndIf
return WavePathList
End
///////////////////////////  END OF GetBrowserSelection Function ///////////////////////////////////////////////////////////


////////////// Sub-Function called from Multiple Functions

Function/WAVE AutExtractWaveFromTable(wRes,Column)
Wave wRes
Variable Column

String name="tempwave"
Variable NPoints=DimSize(wRes,0)

Make/O/N=(NPoints) $name
Wave wname=$name

wname=wRes[p][(Column)]
return wname
End

////////////////////////// End of Extract Column from Table ////////////

/////////////// Function to call Diagnostic Results Function //////////////////////////////////////
/////////// Called from the Macros Menu using the "Check Results of Branch Analysis" Tab
///////// Under Load Data SubMenu
///////// Produces 4 Diagnostic Plots. 3 Show the positions of the ROIs (red) 
/// and corresponding Closest branch nodes (black) in X-Y, X-Z and Y-Z planes
/// The 4th is an histogram of the minimum distances between ROIs and closest node
///// Calls separate function Display_Diagn_Syn for plotting
///// The red arrow that indicates the Current Data Folder in Data Browser needs to point 
///// To the branch you want to analyze
Function Diagnostics()

String name="cell2branch1", Results,SWC
Prompt name, "Name for Branch"
Prompt  Results, "Results Table", popup,Wavelist("*Results",";","")
Prompt  SWC, "SWC File", popup,Wavelist("*swc",";","")
DoPrompt "Input Information", name, Results,SWC

Wave wResults=$Results, wTree=$SWC

Display_Diagn_Syn(wResults,wTree, name)

End
//////////////////////////////////// End call Diagnostic Results Function //////////////////////////////////////
// Displays diagnostic results
Function Display_Diagn_Syn(wResults,wTree, name)
Wave wResults, wTree
String name


// Displays location of nodes and ROIS in 3 separate plots
// X-Y
Display /W=(4.2,39.8,335.4,249.8) wResults[][1] vs wResults[][0] 
AppendtoGraph wResults[][9] vs wResults[][8]
AppendtoGraph wTree[][1] vs wTree[][0]
String tracelist = TraceNameList("", ";", 1), tracename
tracename=StringFromList(1,tracelist)
ModifyGraph mode=3, marker=19, msize=2, rgb($tracename)=(0,0,0)
tracename=StringFromList(2,tracelist)
ModifyGraph mode($tracename)=0,lstyle($tracename)=1,rgb($tracename)=(0,0,0)
Label left "Y co-ordinate (microns)"
Label bottom "X co-ordinate (microns)"
// X-Z
Display /W=(343.8,39.2,721.8,255.8) wResults[][2] vs wResults[][0] 
AppendtoGraph wResults[][10] vs wResults[][8]
AppendtoGraph wTree[][2] vs wTree[][0]
tracelist = TraceNameList("", ";", 1)
tracename=StringFromList(1,tracelist)
ModifyGraph mode=3, marker=19, msize=2, rgb($tracename)=(0,0,0)
tracename=StringFromList(2,tracelist)
ModifyGraph mode($tracename)=0,lstyle($tracename)=1,rgb($tracename)=(0,0,0)
Label left "Z co-ordinate (microns)"
Label bottom "X co-ordinate (microns)"
// Y-Z
Display /W=(6,282.8,338.4,490.4) wResults[][2] vs wResults[][1] 
AppendtoGraph wResults[][10] vs wResults[][9]
AppendtoGraph wTree[][2] vs wTree[][1]
tracelist = TraceNameList("", ";", 1)
tracename=StringFromList(1,tracelist)
ModifyGraph mode=3, marker=19, msize=2, rgb($tracename)=(0,0,0)
tracename=StringFromList(2,tracelist)
ModifyGraph mode($tracename)=0,lstyle($tracename)=1,rgb($tracename)=(0,0,0)
Label left "Z co-ordinate (microns)"
Label bottom "Y co-ordinate (microns)"
// Displays Histogram of Minimum Distances from ROI to node
Make/N=20/O $("Hist"+"MinD"+name)
Wave wHist=$("Hist"+"MinD"+name)

Wave wnodedist=AutExtractWaveFromTable(wResults,7)
Histogram/B=1 wnodedist, wHist
Display /W=(348,283.4,725.4,494) wHist
ModifyGraph mode=5,hbFill=2,rgb=(0,0,0)
Label bottom "Minimum Euclidian 3D Distance to branch (microns)"
KillWaves wnodedist
End
///////////////////////////////	End of Displays diagnostic results ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////// Function to Summarize Results //////////////////////////////////////
/////////// Called from the Macros Menu using the "Summarize Branch Parameters" Tab
///////// Under Analysis SubMenu, uses the SummarizeExc and SummarizeInh functions
/////////// and the Comp_Waves_Dif_Fold_Aut function
///////// Creates Two Tables: One for Excitation (Spine_BrSumm)
////////	and one for Inhibition (Shaft_BrSumm) inside a new Folder called Compile
////// Each Table has a row for each branch
/////// Analyzed columns are: Sum (cumulative Normalized Fluorescence)
////// Avg (Averaged Normalized Fluorescence) Dens (Density synapses per micron)
///// Sum_N (cumulative Normalized Fluorescence Normalize to branch length
///// For Inhibitory table also has Branch Length
//// You can find the column labels switching the horizontal index on the
//// Table Menu to Dimension Labels 
Function Summarize_Sel_Fold()

DFREF RootFolder=GetDataFolderDFR()
String nameTable="Tbl_Inh", nameExc="Summ_Exc", nameInh="Summ_Inh"
String FolderList=GetBrowserSelection_S("")
Variable i=0, L
String dataFolderPathStr, s_i
do
		s_i=num2str(i)
		dataFolderPathStr = StringFromList(i, FolderList, ";")
		if (strlen(dataFolderPathStr) == 0)
			break													// No more data folders.
		endif
		SetDataFolder $dataFolderPathStr 
////////////////////////// Replace here Function to be run in all Folders		
		Wave wInh=$nameTable
		L=wInh[0][2]
		SummarizeExc(L)
		SummarizeInh(L)
		i += 1		
while(1)
String NameCplFolder="Compile", Name_Compile_Exc="Spine_BrSumm", Name_Compile_Inh="Shaft_BrSumm"
Comp_Waves_Dif_Fold_Aut(NameCplFolder,Name_Compile_Exc,nameExc)
Comp_Waves_Dif_Fold_Aut(NameCplFolder,Name_Compile_Inh,nameInh)
SetDataFolder RootFolder
End
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function SummarizeExc(L)
Variable L

String labelSum="Spine_Sum", labelAvg="Spine_Avg", labelDens="Spine_Dens", labelSum_N="Spine_Sum_N"
String name="Summ_Exc",  nameTable="Tbl_Exc"
Make/O/N=(1,4) $name
Wave waveTable=$nameTable

Wave Wexc=AutExtractWaveFromTable(waveTable,1), WRes=$name
WaveStats/Q Wexc
wRes[0][0]=V_Sum
wRes[0][1]=V_avg
wRes[0][2]=V_npnts/L
wRes[0][3]=V_Sum/L
KillWaves Wexc
SetDimLabel 1, 0 , $labelSum,wRes
SetDimLabel 1, 1 , $labelAvg,wRes
SetDimLabel 1, 2 , $labelDens,wRes
SetDimLabel 1, 3 , $labelSum_N,wRes
End
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function SummarizeInh(L)
Variable L

String labelSum="Shaft_Sum", labelAvg="Shaft_Avg", labelDens="Shaft_Dens", labelSum_N="Shaft_Sum_N", labelL="BrLength"
String name="Summ_Inh",  nameTable="Tbl_Inh"
Make/O/N=(1,5) $name
Wave waveTable=$nameTable
Wave WInh_Int=AutExtractWaveFromTable(waveTable,1), WRes=$name
WaveStats/Q WInh_Int
wRes[0][0]=V_Sum
wRes[0][1]=V_avg
wRes[0][2]=V_npnts/L
wRes[0][3]=V_Sum/L
wRes[0][4]=L
KillWaves WInh_Int
SetDimLabel 1, 0 , $labelSum,wRes
SetDimLabel 1, 1 , $labelAvg,wRes
SetDimLabel 1, 2 , $labelDens,wRes
SetDimLabel 1, 3 , $labelSum_N,wRes
SetDimLabel 1, 4 , $labelL,wRes
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////	Compile function Different Folders
Function CompileWavesDifFolders()

String nameFolder, NameWave, NameWaveCompile
Prompt nameFolder, "Name of the Folder to Store Results"
Prompt NameWave, "Name of the Compiled Wave"
Prompt NameWaveCompile, "Name of Wave to Compile"
DoPrompt "Input Variables", nameFolder, NameWave, NameWaveCompile

Comp_Waves_Dif_Fold_aut(nameFolder, NameWave, NameWaveCompile)

End

//////////////////// This is the function called by the Dialogue above
Function Comp_Waves_Dif_Fold_aut(nameFolder, NameWave, NameWaveCompile)
String NameFolder, NameWaveCompile, NameWave



String FolderList=GetBrowserSelection_S("")


String DataFolder=("root:"+NameFolder)

If (DataFolderExists(DataFolder)==0)
	NewDataFolder $DataFolder
else
	SetDataFolder $DataFolder
Endif

String DataFolderPath=DataFolder+":"

Variable i=0
String nextFolderName, s_i, CompileWaveList="", killWaveName
	do
		s_i=num2str(i)
		nextFolderName=StringFromList(i, FolderList)
		if(strlen(nextFolderName)==0)
			break
		endif
		SetDataFolder $nextFolderName
		Wave nextWaveName=$NameWaveCompile
		Duplicate/O nextWaveName, $(NameWaveCompile+s_i)
		Wave dup=$(NameWaveCompile+s_i)
		MoveWave dup, $DataFolderPath
		CompileWaveList+=NameWaveCompile+s_i+";"
		i+=1
	while(1)
SetDataFolder DataFolder
Concatenate/NP=0/O CompileWaveList, $NameWave
i=0
do
	killWaveName=StringFromList(i, CompileWaveList)
	if(strlen(killWaveName)==0)
			break
	endif
	Wave kill=$killWaveName
	KillWaves kill
	i+=1
while(1)

End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function KillWave_Sel_Fold(namewave)
String namewave
// GNeves 05/01/18
// Written as utility
// Removes Named wave in all selected folders
// Assumes Folders to be Analysed are selected in DataBrowser Window
String savedDF= GetDataFolder(1)
/// This Function Must be Present: FolderListFromBrowser()
String FolderList=GetBrowserSelection_S("")

Variable i=0
String dataFolderPathStr
do
		dataFolderPathStr = StringFromList(i, FolderList, ";")
		if (strlen(dataFolderPathStr) == 0)
			break													// No more data folders.
		endif
		SetDataFolder $dataFolderPathStr 
////////////////////////// Replace here Function to be run in all Folders		
		Wave wKill=$namewave
		KillWaves wKill
		i += 1		
while(1)
SetDataFolder root:

End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

