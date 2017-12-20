(* ::Package:: *)

(* ::Section:: *)
(*Begin Package*)


BeginPackage["flowTrack`"];


flowTracks::usage = "flowTracks[{image1,image2..}] takes in a sequence of images and outputs a vector field based on either image correlation or
geometric transformation between sets of points. The window size can be specified as \"windowSize\" -> Integer and
the method as \"method\" -> \"PIV\" or \"PTV\"";


imagePreprocess::usage = "imagePreprocess[image1,image2, opts] takes in two images and a list of preprocessing options to process images"


PIV::usage = "PIV[image1,image2,windowsize,method] takes in an image pair, a window size and a method to perform a simple PIV"


PTV::usage = "PTV[image1,image2,windowsize] yields a crude PTV strategy"


(* ::Section:: *)
(*Main Functions*)


Begin["`Private`"];


Options[imagePreprocess]={"clipped"-> False,"histogramEqualize"-> False,"highP"-> False,"wiener"-> False,"wienerSize"-> 3, "highPSize"-> 15}
imagePreprocess[image1_, image2_,OptionsPattern[]]:=
Block[{img ,upperlimit},
img = ColorConvert[#,"Grayscale"]&/@{image1,image2};
If[OptionValue["clipped"], 
upperlimit = Map[Median[#]+2 StandardDeviation[#]&,img];
img = MapThread[ReplacePixelValue[#1,
SparseArray[ImageData@Clip[#1-#2]]["NonzeroPositions"]->#2]&,
{img,upperlimit}]
];
If[OptionValue["histogramEqualize"], img = HistogramTransform/@img];
If[OptionValue["highP"], img = Map[(# - GaussianFilter[#,ConstantArray[OptionValue["highPSize"],2]])/Max[#]*255 &, img]];
If[OptionValue["wiener"], img = WienerFilter[#,OptionValue["wienerSize"]]&/@img];
img
];


PIV[image1_?ImageQ,image2_?ImageQ,win_Integer, pivmethod_]:=Module[{windowsize=win,imgDim=ImageDimensions[image1],
img1NoBorder,interrogateWin,searchWins,dim, midPtsImg1,correlationPts,h,f,width,img1,img2,imgdata,height},

{img1,img2} = imagePreprocess[image1,image2]; (* image preprocessing *)

{width,height}=imgDim;
img1NoBorder=ImageCrop[img1,imgDim-(2*windowsize)]; (*removing border from the first image*)
interrogateWin=img1NoBorder~ImagePartition~windowsize;(*splitting the image (without borders) into small windows. 
correlations for each  of these windows will be found in the second image*)
searchWins=ImagePartition[img2,3*windowsize,{windowsize,windowsize}];
(*breaking the second image into search Windows. Each search window has a dimension three times that of the entry in interrogateWin*)
dim=Dimensions[searchWins];
{f,h}= ImageDimensions[img1NoBorder];(*finding the midpoints of the interrogateWin*)
midPtsImg1= Transpose[Table[Abs[{i windowsize+windowsize/2,j(windowsize)+windowsize/2}-{0,height}],{i,Last@dim},{j,First@dim}]]~Flatten~1;
(*finding correlation points between interrogateWin (splitwindows of the first image) and searchWins (search windows) of the
second image*)
correlationPts=Table[
imgdata=ImageData@interrogateWin[[i+1,j+1]];
If[And@@Map[First[imgdata]==#&,imgdata,{1}],
Indeterminate,
Composition[
{#[[2]],height-#[[1]]}&,
FirstPosition[0],
MorphologicalComponents]@ImagePad[ImageAdjust@ImageCorrelate[searchWins[[i+1,j+1]],interrogateWin[[i+1,j+1]],
NormalizedSquaredEuclideanDistance,PerformanceGoal->"Quality"],{{j windowsize ,f-windowsize(j+1)},{h-windowsize(i+1),i windowsize}},
White]],{i,0,First[dim]-1},{j,0,Last[dim]-1}];
correlationPts=Cases[correlationPts,{__Integer}|Indeterminate,Infinity];
{midPtsImg1,correlationPts}=DeleteCases[Thread[{midPtsImg1,correlationPts}],{_,Indeterminate}]\[Transpose];
(*sow the graphics window with the vectors*)
Sow@Graphics[{Arrowheads[0.01],MapThread[Arrow[{##}]&,{midPtsImg1,correlationPts}]},Frame->True];
image2
];


PTV[image1_?ImageQ,image2_?ImageQ,win_Integer]:=Module[{imgDim=ImageDimensions[image1],imgCorrD,img,geomtranF,rpts},
imgCorrD=First@imgDim;
img=Map[Binarize]@(ImageCrop[#,imgDim-(2*win)]&/@{image1,image2});(*binarize the images*)
rpts=PixelValuePositions[ImagePad[First@img,win],1];(*detecting the particles in the first image*)
geomtranF=Last[FindGeometricTransform@@img];(*finding the geometric transform*)
(*sow the vector-field*)
Sow@Graphics[{Red,Arrowheads[.01], Arrow/@Transpose[{Map[Abs[#-{0,imgCorrD}]&,rpts],Map[Abs[#-{0,imgCorrD}]&,geomtranF[rpts]]}]},Frame->True];
image2
]


Options[flowTracks]={"method"->"PIV","windowSize"->32,"pivmethod"-> NormalizedSquaredEuclideanDistance};
flowTracks[p:{images__?ImageQ},OptionsPattern[]]:=Module[{opt=OptionValue["method"], win=OptionValue["windowSize"],pivmethod = OptionValue["pivmethod"]},
Switch[opt,"PIV",
(Last@Reap@FoldList[PIV[#1,#2,win,pivmethod]&,First[p],Rest[p]])~Flatten~1,
"PTV",
(Last@Reap@FoldList[PTV[#1,#2,win]&,First[p],Rest[p]])~Flatten~1
]
];


(* ::Section:: *)
(*End Package*)


End[];


EndPackage[];
