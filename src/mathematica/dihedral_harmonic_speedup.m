(* Content-type: application/vnd.wolfram.mathematica *)

(*** Wolfram Notebook File ***)
(* http://www.wolfram.com/nb *)

(* CreatedBy='Mathematica 13.3' *)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[       158,          7]
NotebookDataLength[     41843,        923]
NotebookOptionsPosition[     35791,        835]
NotebookOutlinePosition[     36229,        853]
CellTagsIndexPosition[     36186,        850]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Initialization code", "Title",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{
  3.402914001047615*^9, {3.952051117940569*^9, 3.9520511268304253`*^9}, 
   3.952207602521706*^9, 
   3.952350500210045*^9},ExpressionUUID->"59ad4f77-f93b-4129-8c5d-\
5c71330c705c"],

Cell[BoxData[
 RowBox[{"$RecursionLimit", " ", "=", " ", 
  RowBox[{"Floor", "[", 
   RowBox[{"(", 
    RowBox[{"2", "^", "12"}], ")"}], "]"}]}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.910222496749535*^9, 3.910222502985106*^9}, 
   3.952207602521844*^9, {3.9523450606915617`*^9, 3.9523450627582493`*^9}, {
   3.952345805604395*^9, 3.9523458066264877`*^9}, {3.952345943503106*^9, 
   3.952345944778145*^9}, {3.952346014592503*^9, 3.952346039047865*^9}, 
   3.952347762306674*^9, 3.9523505002111483`*^9},
 CellLabel->"In[1]:=",ExpressionUUID->"2302b7fc-c3e0-473a-9944-bd8907480542"],

Cell[BoxData[
 RowBox[{"Directory", "[", "]"}]], "Code",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.402675593015943*^9, 3.40267560234683*^9}, 
   3.402675719919037*^9, 3.402914001047728*^9, 3.952207602521925*^9, 
   3.9523505002119417`*^9},
 CellLabel->"In[2]:=",ExpressionUUID->"571b4188-8844-47fd-a5e9-9e8bc63da377"],

Cell[BoxData[
 RowBox[{"SetDirectory", "[", 
  "\"\</Users/meister/Development/cando/extensions/cando/include/cando/chem/\
energy_functions/\>\"", "]"}]], "Code",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.402675724744155*^9, 3.402675744579932*^9}, {
   3.4026757812426653`*^9, 3.402675795163514*^9}, 3.4029140010477962`*^9, {
   3.404776927640129*^9, 3.4047769282020483`*^9}, {3.485009516098822*^9, 
   3.4850095173434477`*^9}, {3.619556879758088*^9, 3.6195568829874687`*^9}, {
   3.7026648712565193`*^9, 3.702664872815674*^9}, {3.7809265423187447`*^9, 
   3.780926543584301*^9}, {3.780944859091173*^9, 3.780944904216384*^9}, {
   3.8946458653781567`*^9, 3.894645866645886*^9}, 3.9522076025221043`*^9, 
   3.952350500212072*^9},
 CellLabel->"In[3]:=",ExpressionUUID->"13260320-554d-4f6a-97fe-a104adcbae00"],

Cell["Setup to generate code that is embedded within CANDO", "Subtitle",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{
  3.402914001047846*^9, {3.494842902549116*^9, 3.494842937561069*^9}, 
   3.894645879760284*^9, {3.951978343620945*^9, 3.951978359627625*^9}, {
   3.951979633876155*^9, 3.9519796348630543`*^9}, 
   3.952350500212597*^9},ExpressionUUID->"f3c71499-dcd8-4a60-8091-\
9938479e0665"],

Cell[BoxData[
 RowBox[{"NotebookDirectory", "[", "]"}]], "Code",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.951979559481789*^9, 3.951979570445331*^9}, 
   3.952350500212707*^9},
 CellLabel->"In[4]:=",ExpressionUUID->"29cf2b5b-9ced-479c-85ae-c0a9be6c68e4"],

Cell[BoxData[
 RowBox[{"$Path", "=", 
  RowBox[{"Append", "[", 
   RowBox[{"$Path", ",", 
    RowBox[{"NotebookDirectory", "[", "]"}]}], "]"}]}]], "Code",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.9519795768887367`*^9, 3.951979587329652*^9}, 
   3.952350500212812*^9},
 CellLabel->"In[5]:=",ExpressionUUID->"8c98e9fb-8d7f-4365-9e46-d53e1da154f8"],

Cell[BoxData[
 RowBox[{"Needs", "[", "\"\<foo`\>\"", "]"}]], "Code",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.4055916304006453`*^9, 3.405591641125656*^9}, 
   3.9523505002129183`*^9, {3.952944081092002*^9, 3.95294408160835*^9}, 
   3.952945905737265*^9, {3.952947930867947*^9, 3.95294793203485*^9}, {
   3.95294812546375*^9, 3.9529481256572857`*^9}},
 CellLabel->"In[6]:=",ExpressionUUID->"420f0456-2a09-4d78-bd8d-6969828b5b21"]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[{
 "DihedralHarmonic term using ",
 StyleBox["Mathematica",
  FontSlant->"Italic"],
 " generated gradient/Hessian expansion"
}], "Chapter",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.494844035232679*^9, 3.494844038080626*^9}, 
   3.952260616260887*^9, 3.952350500213327*^9, {3.9525115463144693`*^9, 
   3.952511548817813*^9}},ExpressionUUID->"1f1aa4cc-2a9c-461d-93ec-\
97c2242c265e"],

Cell[BoxData[
 RowBox[{
  RowBox[{"Clear", "[", 
   RowBox[{
   "x", ",", "y", ",", "z", ",", "x1", ",", "y1", ",", "z1", ",", "x2", ",", 
    "y2", ",", "z2", ",", "x3", ",", "y3", ",", "z3", ",", "x4", ",", "y4", 
    ",", "z4", ",", "pa", ",", "pb", ",", "pc", ",", "pd", ",", "kdh", ",", 
    " ", "phi0"}], "]"}], ";"}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.402652780531734*^9, 3.402652792729147*^9}, {
   3.405593978391768*^9, 3.405593979390263*^9}, {3.9522074423482924`*^9, 
   3.952207487258544*^9}, {3.9523331896757097`*^9, 3.952333196002604*^9}, {
   3.95233346309575*^9, 3.9523334642368393`*^9}, 3.952350500214826*^9},
 CellLabel->"In[7]:=",ExpressionUUID->"eda73ac1-ca89-4acb-93f8-a2bcaef0ea63"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"pa", " ", "=", " ", 
   RowBox[{"{", " ", 
    RowBox[{"x1", ",", " ", "y1", ",", " ", "z1"}], "}"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"pb", " ", "=", " ", 
   RowBox[{"{", " ", 
    RowBox[{"x2", ",", " ", "y2", ",", " ", "z2"}], " ", "}"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"pc", " ", "=", " ", 
   RowBox[{"{", " ", 
    RowBox[{"x3", ",", " ", "y3", ",", " ", "z3"}], " ", "}"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"pd", " ", "=", " ", 
   RowBox[{"{", " ", 
    RowBox[{"x4", ",", " ", "y4", ",", " ", "z4"}], " ", "}"}]}], 
  ";"}]}], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.40265249626128*^9, 3.4026525067255383`*^9}, {
   3.402652543649145*^9, 3.4026526050967493`*^9}, {3.9522070732514353`*^9, 
   3.952207081595791*^9}, {3.9522075124882*^9, 3.952207518226007*^9}, 
   3.9523505002166557`*^9},
 CellLabel->"In[8]:=",ExpressionUUID->"ab8ce859-bbb6-4364-b2a0-760d04215a00"],

Cell[BoxData[
 RowBox[{
  RowBox[{"deg", "[", "r_", "]"}], " ", ":=", " ", 
  RowBox[{"r", "/", "0.0174533"}]}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.952112299136571*^9, 3.9521123134850273`*^9}, 
   3.952350500216721*^9},
 CellLabel->"In[12]:=",ExpressionUUID->"f4046a73-8241-4cc3-92b5-b0168510db9b"],

Cell[BoxData[
 RowBox[{
  RowBox[{"rad", "[", "d_", "]"}], " ", ":=", " ", 
  RowBox[{"d", "*", "0.0174533"}]}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.95211789101343*^9, 3.9521179011880217`*^9}, 
   3.952350500216776*^9},
 CellLabel->"In[13]:=",ExpressionUUID->"eda20186-450b-48df-880f-677c611f8848"],

Cell[BoxData[
 RowBox[{
  RowBox[{"dihedralAngle", "[", 
   RowBox[{"{", 
    RowBox[{"x1_", ",", "x2_", ",", "x3_", ",", "x4_"}], "}"}], "]"}], " ", ":=",
   "\[IndentingNewLine]", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
     "d1", ",", "d2", ",", "d3", ",", "b1", ",", "b2", ",", "b3", ",", "n1", 
      ",", "n2", ",", "m1", ",", "x", ",", "y", ",", "theta"}], "}"}], ",", 
    "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"d1", " ", "=", " ", 
      RowBox[{"(", 
       RowBox[{"x2", "-", "x1"}], ")"}]}], ";", "\[IndentingNewLine]", 
     RowBox[{"d2", " ", "=", " ", 
      RowBox[{"(", 
       RowBox[{"x3", "-", "x2"}], ")"}]}], ";", "\[IndentingNewLine]", 
     RowBox[{"d3", " ", "=", " ", 
      RowBox[{"(", 
       RowBox[{"x4", "-", "x3"}], ")"}]}], ";", "\[IndentingNewLine]", 
     RowBox[{"b1", "=", 
      RowBox[{"d1", "/", 
       RowBox[{"Sqrt", "[", 
        RowBox[{"Dot", "[", 
         RowBox[{"d1", ",", "d1"}], "]"}], "]"}]}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"b2", "=", 
      RowBox[{"d2", "/", 
       RowBox[{"Sqrt", "[", 
        RowBox[{"Dot", "[", 
         RowBox[{"d2", ",", "d2"}], "]"}], "]"}]}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"b3", "=", 
      RowBox[{"d3", "/", 
       RowBox[{"Sqrt", "[", 
        RowBox[{"Dot", "[", 
         RowBox[{"d3", ",", "d3"}], "]"}], "]"}]}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"n1", "=", 
      RowBox[{"Cross", "[", 
       RowBox[{"b1", ",", "b2"}], "]"}]}], ";", "\[IndentingNewLine]", 
     RowBox[{"n2", "=", 
      RowBox[{"Cross", "[", 
       RowBox[{"b2", ",", "b3"}], "]"}]}], ";", "\[IndentingNewLine]", 
     RowBox[{"m1", "=", 
      RowBox[{"Cross", "[", 
       RowBox[{"n1", ",", "b2"}], "]"}]}], ";", "\[IndentingNewLine]", 
     RowBox[{"x", "=", 
      RowBox[{"n1", ".", "n2"}]}], ";", "\[IndentingNewLine]", 
     RowBox[{"y", "=", 
      RowBox[{"m1", ".", "n2"}]}], ";", "\[IndentingNewLine]", 
     RowBox[{"theta", "=", 
      RowBox[{"ArcTan", "[", 
       RowBox[{"x", ",", 
        RowBox[{"-", "y"}]}], "]"}]}], ";", "\[IndentingNewLine]", 
     "theta"}]}], "\[IndentingNewLine]", "]"}]}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.952115964310054*^9, 3.952116038005422*^9}, {
   3.95211608996492*^9, 3.952116098316346*^9}, {3.952116174490715*^9, 
   3.952116176524521*^9}, {3.952116209027203*^9, 3.952116214909787*^9}, 
   3.952116255016013*^9, 3.952116311243758*^9, 3.95211663078019*^9, 
   3.952117050351493*^9, 3.952350500215158*^9, {3.952352634046039*^9, 
   3.952352738955731*^9}, {3.9523599419921503`*^9, 3.952359950737565*^9}, {
   3.9523603299297*^9, 3.9523603335353317`*^9}, 3.952365230057803*^9},
 CellLabel->"In[14]:=",ExpressionUUID->"9d21da8e-a00b-4bbf-819a-d0b689cd5914"],

Cell[BoxData[
 RowBox[{
  RowBox[{"phiFn", " ", "=", " ", 
   RowBox[{"dihedralAngle", "[", 
    RowBox[{"{", 
     RowBox[{"pa", ",", "pb", ",", "pc", ",", "pd"}], "}"}], "]"}]}], 
  ";"}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.952369378010898*^9, 3.95236938983106*^9}, 
   3.952511853098412*^9, 3.952511978521289*^9},
 CellLabel->"In[15]:=",ExpressionUUID->"84b668a5-368f-441a-98ba-2cff9272ae57"],

Cell[BoxData[
 RowBox[{
  RowBox[{"deltaPhiModFn", " ", "=", " ", 
   RowBox[{
    RowBox[{"Mod", "[", 
     RowBox[{
      RowBox[{"phiFn", "-", "phi0", "+", "Pi"}], ",", 
      RowBox[{"2.0", "*", "Pi"}]}], "]"}], "-", "Pi"}]}], ";"}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.952364231952064*^9, 3.9523642908052263`*^9}, {
  3.9523658996415577`*^9, 3.95236592176656*^9}, {3.9523693951854753`*^9, 
  3.952369396076029*^9}, {3.95251198279893*^9, 3.952511993614575*^9}},
 CellLabel->"In[16]:=",ExpressionUUID->"115dfce1-7a62-4a15-8b38-98017ab13dfc"],

Cell[BoxData[
 RowBox[{
  RowBox[{"dihedralHarmonicInputs", " ", "=", " ", 
   RowBox[{"Flatten", "[", 
    RowBox[{"{", 
     RowBox[{"pa", ",", "pb", ",", "pc", ",", "pd", ",", 
      RowBox[{"{", 
       RowBox[{"phi0", ",", "kdh"}], "}"}]}], "}"}], "]"}]}], ";"}]], "Input",\

 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.4055966209501877`*^9, 3.405596670762001*^9}, {
   3.405596704315874*^9, 3.405596706769898*^9}, {3.410178228193839*^9, 
   3.410178228313593*^9}, {3.9522070926667356`*^9, 3.952207123863529*^9}, {
   3.952207262876726*^9, 3.952207286428831*^9}, {3.95220752424325*^9, 
   3.952207532509082*^9}, 3.952333205148869*^9, 3.952350500216949*^9, 
   3.952364795153647*^9},
 CellLabel->"In[17]:=",ExpressionUUID->"f7032a5f-8737-4c9b-9eca-017969e08d46"],

Cell[BoxData[
 RowBox[{
  RowBox[{"dihedralHarmonicEnergyFn", " ", "=", " ", 
   RowBox[{
    FractionBox["1", "2"], "kdh", 
    SuperscriptBox[
     RowBox[{"(", "deltaPhiModFn", ")"}], "2"]}]}], ";"}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.952052587078174*^9, 3.952052631968499*^9}, {
   3.952052745649139*^9, 3.952052757535961*^9}, {3.952061993482945*^9, 
   3.952062057594069*^9}, {3.952062120767469*^9, 3.952062124086547*^9}, {
   3.952094250875205*^9, 3.952094264376678*^9}, {3.952333228583371*^9, 
   3.952333245708335*^9}, 3.952350500218419*^9, {3.952360096657234*^9, 
   3.952360115334506*^9}, {3.952361690302788*^9, 3.9523617245340014`*^9}, {
   3.952364408072691*^9, 3.952364412088141*^9}, {3.95236444596378*^9, 
   3.9523644559496107`*^9}, 3.9523645011948347`*^9, {3.952364585521674*^9, 
   3.952364586266634*^9}, 3.952364764759695*^9, 3.952364807765966*^9, {
   3.952365934346266*^9, 3.9523659624388533`*^9}, {3.95236635926016*^9, 
   3.952366371136812*^9}, {3.9525120289973288`*^9, 3.952512046567223*^9}, {
   3.9525121736339617`*^9, 3.952512180637619*^9}, {3.952512301041863*^9, 
   3.9525123145494947`*^9}},
 CellLabel->"In[18]:=",ExpressionUUID->"c3d70802-8e89-4cfd-aeb7-755b7ecd4da5"],

Cell[BoxData[
 RowBox[{"dihedralHarmonicEnergyFromDeltaPhiVarFn", " ", "=", " ", 
  RowBox[{"kdh", 
   RowBox[{
    RowBox[{"(", "DeltaPhiVar", ")"}], "^", "2"}]}]}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.9523654052424297`*^9, 3.952365485642376*^9}, 
   3.9523663641527157`*^9, {3.9525122731675777`*^9, 3.952512285893382*^9}},
 CellLabel->"In[19]:=",ExpressionUUID->"d4b06aff-2869-4380-8cb9-fe70fe3b34b8"],

Cell[BoxData[
 RowBox[{
  RowBox[{"parms", "[", 
   RowBox[{"phi_", ",", "xphi0_", ",", "xkdh_"}], "]"}], " ", ":=", 
  RowBox[{"{", "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"kdh", "->", "xkdh"}], ",", 
    RowBox[{"phi0", "->", "xphi0"}], ",", "\[IndentingNewLine]", 
    RowBox[{"x1", "->", "0"}], ",", 
    RowBox[{"y1", "->", "1"}], ",", 
    RowBox[{"z1", "->", "0"}], ",", "\[IndentingNewLine]", 
    RowBox[{"x2", "->", "0"}], ",", 
    RowBox[{"y2", "->", "0"}], ",", 
    RowBox[{"z2", "->", "0"}], ",", "\[IndentingNewLine]", 
    RowBox[{"x3", "->", "1"}], ",", 
    RowBox[{"y3", "->", "0"}], ",", 
    RowBox[{"z3", "->", "0"}], ",", "\[IndentingNewLine]", 
    RowBox[{"x4", "->", "1"}], ",", 
    RowBox[{"y4", "->", 
     RowBox[{"Cos", "[", "phi", "]"}]}], ",", 
    RowBox[{"z4", "->", 
     RowBox[{"Sin", "[", "phi", "]"}]}]}], "}"}]}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.952364979921091*^9, 3.95236510763378*^9}, 
   3.9523651441804543`*^9, 3.952365247996396*^9},
 CellLabel->"In[20]:=",ExpressionUUID->"8b9e4bcf-6850-4d82-bade-26d978769d42"],

Cell[BoxData[
 RowBox[{"Plot", "[", 
  RowBox[{
   RowBox[{"dihedralHarmonicEnergyFn", "/.", 
    RowBox[{"parms", "[", 
     RowBox[{"phi", ",", 
      RowBox[{"0.5", "*", "Pi"}], ",", "100"}], "]"}]}], ",", 
   RowBox[{"{", 
    RowBox[{"phi", ",", 
     RowBox[{"-", "Pi"}], ",", "Pi"}], "}"}], ",", 
   RowBox[{"PlotRange", "->", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{"-", "Pi"}], ",", "Pi"}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"0", ",", "10"}], "}"}]}], "}"}]}]}], "]"}]], "Input",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.},
 CellChangeTimes->{{3.952366140597803*^9, 3.952366197957847*^9}, {
  3.95251214149334*^9, 3.952512142818516*^9}, {3.952512618586171*^9, 
  3.952512640208448*^9}},
 CellLabel->"In[21]:=",ExpressionUUID->"34e7ce02-3852-496d-9c58-948a1952826d"],

Cell[BoxData[
 RowBox[{"D", "[", 
  RowBox[{"dihedralHarmonicEnergyFromDeltaPhiVarFn", ",", "DeltaPhiVar"}], 
  "]"}]], "Input",
 CellChangeTimes->{
  3.9523653648522797`*^9, {3.9523654643077097`*^9, 3.952365502657693*^9}, 
   3.952512417327572*^9},
 CellLabel->"In[22]:=",ExpressionUUID->"266cf776-52f1-4a43-b178-bdf8459ae224"],

Cell[BoxData[
 RowBox[{"names", " ", "=", " ", 
  RowBox[{"Flatten", "[", 
   RowBox[{"{", 
    RowBox[{"pa", ",", "pb", ",", "pc", ",", "pd"}], "}"}], "]"}]}]], "Input",\

 CellChangeTimes->{{3.406839124224272*^9, 3.406839125529244*^9}, {
  3.4068391856032763`*^9, 3.406839188261937*^9}, {3.952366251486764*^9, 
  3.952366254037086*^9}},
 CellLabel->"In[23]:=",ExpressionUUID->"590d8c9e-4792-413c-9810-3bb0c7f7e9a1"],

Cell["Energy", "Subsection",
 CellChangeTimes->{{3.4017146731336727`*^9, 
  3.4017146744551163`*^9}},ExpressionUUID->"6fb2043b-fed8-4159-9046-\
f2a26ced5774"],

Cell[BoxData[
 RowBox[{
  RowBox[{"dihedralHarmonicVarNames", " ", "=", " ", 
   RowBox[{"{", " ", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"x1", ",", "\"\<x\>\"", ",", "I1", ",", "0"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"y1", ",", "\"\<y\>\"", ",", "I1", ",", "1"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"z1", ",", "\"\<z\>\"", ",", "I1", ",", "2"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"x2", ",", "\"\<x\>\"", ",", "I2", ",", "0"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"y2", ",", "\"\<y\>\"", ",", "I2", ",", "1"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"z2", ",", "\"\<z\>\"", ",", "I2", ",", "2"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"x3", ",", "\"\<x\>\"", ",", "I3", ",", "0"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"y3", ",", "\"\<y\>\"", ",", "I3", ",", "1"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"z3", ",", "\"\<z\>\"", ",", "I3", ",", "2"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"x4", ",", "\"\<x\>\"", ",", "I4", ",", "0"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"y4", ",", "\"\<y\>\"", ",", "I4", ",", "1"}], "}"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"z4", ",", "\"\<z\>\"", ",", "I4", ",", "2"}], "}"}]}], 
    "\[IndentingNewLine]", "}"}]}], ";"}]], "Input",
 CellChangeTimes->{{3.4085450415029783`*^9, 3.408545189001713*^9}, {
   3.40854564425465*^9, 3.40854566509382*^9}, {3.408545841159666*^9, 
   3.40854589498903*^9}, {3.408562995604747*^9, 3.408563013307331*^9}, {
   3.95233495193824*^9, 3.952334980758029*^9}, 3.952365629726313*^9, {
   3.952366260800379*^9, 3.952366273353196*^9}},
 CellLabel->"In[24]:=",ExpressionUUID->"42d21f14-444c-49d2-9590-1099b09c7d89"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"dihedralHarmonicSetupRules", " ", "=", " ", 
   RowBox[{"{", "}"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicSetupRules", ",", 
    RowBox[{
    "CCode", "[", "\"\<DIHEDRAL_HARMONIC_SET_PARAMETER(kdh);\>\"", "]"}]}], 
   "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicSetupRules", ",", 
    RowBox[{
    "CCode", "[", "\"\<DIHEDRAL_HARMONIC_SET_PARAMETER(phi0);\>\"", "]"}]}], 
   "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicSetupRules", ",", 
    RowBox[{
    "CCode", "[", "\"\<DIHEDRAL_HARMONIC_SET_PARAMETER(I1);\>\"", "]"}]}], 
   "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicSetupRules", ",", 
    RowBox[{
    "CCode", "[", "\"\<DIHEDRAL_HARMONIC_SET_PARAMETER(I2);\>\"", "]"}]}], 
   "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicSetupRules", ",", 
    RowBox[{
    "CCode", "[", "\"\<DIHEDRAL_HARMONIC_SET_PARAMETER(I3);\>\"", "]"}]}], 
   "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicSetupRules", ",", 
    RowBox[{
    "CCode", "[", "\"\<DIHEDRAL_HARMONIC_SET_PARAMETER(I4);\>\"", "]"}]}], 
   "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicSetupRules", ",", 
    RowBox[{
    "CCode", "[", "\"\<DIHEDRAL_HARMONIC_APPLY_ATOM_MASK(I1,I2,I3,I4);\>\"", 
     "]"}]}], "]"}], ";"}]}], "Input",
 CellChangeTimes->{{3.410177124188498*^9, 3.410177277864727*^9}, {
  3.4101777777433367`*^9, 3.410177805652183*^9}, {3.410178179796535*^9, 
  3.410178203278604*^9}, {3.4101782360285482`*^9, 3.410178236128345*^9}, {
  3.910226587560244*^9, 3.910226605409601*^9}, {3.952365630850531*^9, 
  3.952365654280559*^9}, {3.952365685029888*^9, 3.952365705308198*^9}, {
  3.9523657796346283`*^9, 3.9523658208752813`*^9}},
 CellLabel->"In[25]:=",ExpressionUUID->"a18728c5-d7d6-43e3-8ba0-e33e1aeffaa8"],

Cell[BoxData[
 RowBox[{
  RowBox[{"For", "[", 
   RowBox[{
    RowBox[{"i", "=", "1"}], ",", 
    RowBox[{"i", "\[LessEqual]", 
     RowBox[{"Length", "[", "dihedralHarmonicVarNames", "]"}]}], ",", 
    RowBox[{"i", "++"}], ",", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"str", " ", "=", " ", 
      RowBox[{"\"\<DIHEDRAL_HARMONIC_SET_POSITION(\>\"", "<>", 
       RowBox[{"ToString", "[", 
        RowBox[{
         RowBox[{"dihedralHarmonicVarNames", "[", 
          RowBox[{"[", "i", "]"}], "]"}], "[", 
         RowBox[{"[", "1", "]"}], "]"}], "]"}], "<>", "\"\<,\>\"", "<>", 
       RowBox[{"ToString", "[", 
        RowBox[{
         RowBox[{"dihedralHarmonicVarNames", "[", 
          RowBox[{"[", "i", "]"}], "]"}], "[", 
         RowBox[{"[", "3", "]"}], "]"}], "]"}], "<>", "\"\<,\>\"", "<>", 
       RowBox[{"ToString", "[", 
        RowBox[{
         RowBox[{"dihedralHarmonicVarNames", "[", 
          RowBox[{"[", "i", "]"}], "]"}], "[", 
         RowBox[{"[", "4", "]"}], "]"}], "]"}], "<>", "\"\<);\>\""}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"PrintTemporary", "[", "str", "]"}], ";", "\[IndentingNewLine]", 
     
     RowBox[{"AppendTo", "[", 
      RowBox[{"dihedralHarmonicSetupRules", ",", 
       RowBox[{"CCode", "[", "str", "]"}]}], "]"}], ";"}]}], 
   "\[IndentingNewLine]", "]"}], ";"}]], "Input",
 CellChangeTimes->{{3.410177280324893*^9, 3.410177373744402*^9}, {
   3.410177477477*^9, 3.410177533673448*^9}, {3.410177572074147*^9, 
   3.4101775752355337`*^9}, {3.410177780093665*^9, 3.410177806677978*^9}, {
   3.9523656553573217`*^9, 3.952365658838181*^9}, 3.9523657059805098`*^9, {
   3.952366708899272*^9, 3.952366711322466*^9}},
 CellLabel->"In[33]:=",ExpressionUUID->"db4a22fd-feda-40c4-942d-77cc9d3ce06a"],

Cell[BoxData[
 RowBox[{"dihedralHarmonicSetupRules", "//", "MatrixForm"}]], "Input",
 CellChangeTimes->{{3.410177392336293*^9, 3.410177395110338*^9}, {
   3.41017756267308*^9, 3.4101775644454517`*^9}, 3.410177784928906*^9, 
   3.952365660040989*^9},
 CellLabel->"In[34]:=",ExpressionUUID->"9888ee28-cae6-446b-a9ed-7bb0a5b91811"],

Cell[BoxData[
 RowBox[{
  RowBox[{"dihedralHarmonicOutputs", " ", "=", " ", 
   RowBox[{"{", "}"}]}], ";"}]], "Input",
 CellChangeTimes->{{3.408587821227968*^9, 3.4085878252235737`*^9}, 
   3.9523656608508797`*^9},
 CellLabel->"In[35]:=",ExpressionUUID->"c1775d99-a2a5-4151-95e4-b6accbfb9c97"],

Cell[BoxData[
 RowBox[{"dihedralHarmonicHessian", " ", "=", " ", 
  RowBox[{"Table", "[", 
   RowBox[{
    RowBox[{"Table", "[", 
     RowBox[{"0", ",", 
      RowBox[{"{", "12", "}"}]}], "]"}], ",", 
    RowBox[{"{", "12", "}"}]}], "]"}]}]], "Input",
 CellChangeTimes->{{3.40933335742642*^9, 3.40933336900301*^9}, 
   3.952365661419042*^9, {3.9523658293207197`*^9, 3.952365830771131*^9}},
 CellLabel->"In[36]:=",ExpressionUUID->"d45175be-4779-4d6a-ba17-88e005b2e201"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"dihedralHarmonicEnergyRules", " ", "=", " ", 
   RowBox[{"{", "}"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicEnergyRules", ",", 
    RowBox[{"Assign", "[", 
     RowBox[{"DeltaPhiVar", ",", "deltaPhiModFn"}], "]"}]}], "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicEnergyRules", ",", 
    RowBox[{"Assign", "[", 
     RowBox[{"Energy", ",", "dihedralHarmonicEnergyFromDeltaPhiVarFn"}], 
     "]"}]}], "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicEnergyRules", ",", 
    RowBox[{"EnergyAccumulate", "[", 
     RowBox[{"\"\<DIHEDRAL_HARMONIC\>\"", ",", "Energy"}], "]"}]}], "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicOutputs", ",", "DeltaPhiVar"}], "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendTo", "[", 
   RowBox[{"dihedralHarmonicOutputs", ",", "Energy"}], "]"}], ";"}]}], "Input",\

 CellChangeTimes->{{3.401714154505054*^9, 3.40171415874859*^9}, {
   3.401714286121945*^9, 3.401714287127378*^9}, {3.401714414109407*^9, 
   3.4017144145268173`*^9}, {3.402652075438241*^9, 3.402652079400985*^9}, {
   3.402660211700779*^9, 3.402660211962166*^9}, 3.402660717868671*^9, {
   3.403183975900215*^9, 3.403183982112797*^9}, {3.405591741230144*^9, 
   3.405591745357368*^9}, {3.405605212352811*^9, 3.405605255675804*^9}, {
   3.405606432112876*^9, 3.405606432293248*^9}, {3.408552349903845*^9, 
   3.408552365804942*^9}, {3.408553616504222*^9, 3.408553632331004*^9}, {
   3.40855893637739*^9, 3.408558954320676*^9}, {3.408563249888267*^9, 
   3.408563318638526*^9}, {3.4085635055620623`*^9, 3.408563506089058*^9}, {
   3.409106851769322*^9, 3.409106882056981*^9}, {3.409333352642502*^9, 
   3.40933335434091*^9}, {3.4101778333157253`*^9, 3.410177835865591*^9}, {
   3.410177870714059*^9, 3.410177885725202*^9}, {3.410366745442485*^9, 
   3.4103668611031027`*^9}, {3.410387627352058*^9, 3.410387662996006*^9}, {
   3.410388393326111*^9, 3.4103883934847193`*^9}, {3.4103920130578737`*^9, 
   3.410392015575086*^9}, {3.410392087452957*^9, 3.4103920915966787`*^9}, {
   3.410392140741385*^9, 3.41039215558921*^9}, 3.410743057646208*^9, {
   3.411174874722517*^9, 3.411174878860882*^9}, {3.411174974239976*^9, 
   3.4111749916562967`*^9}, {3.41119932577212*^9, 3.4111993265205603`*^9}, {
   3.9102765561456757`*^9, 3.910276558525859*^9}, {3.9523656629081097`*^9, 
   3.9523656669440813`*^9}, 3.9523657339789248`*^9, {3.952365765653974*^9, 
   3.952365767244935*^9}, 3.952365839304942*^9, {3.9523663004005632`*^9, 
   3.952366342358136*^9}, 3.952366413351788*^9, 3.952366511638007*^9, 
   3.952366554255515*^9, {3.952368297620332*^9, 3.952368310533153*^9}},
 CellLabel->"In[37]:=",ExpressionUUID->"fa831a18-e7d4-42d3-86d1-85c6f7869e16"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"dihedralHarmonicForceHessianRules", " ", "=", " ", 
   RowBox[{"{", "}"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"AppendChainGradientForceAndHessian", "[", "\[IndentingNewLine]", 
   RowBox[{
   "\"\<DIHEDRAL_HARMONIC\>\"", ",", "\[IndentingNewLine]", 
    "dihedralHarmonicForceHessianRules", ",", "\[IndentingNewLine]", 
    "dihedralHarmonicOutputs", ",", "\[IndentingNewLine]", 
    "dihedralHarmonicHessian", ",", "\[IndentingNewLine]", 
    "dihedralHarmonicEnergyFromDeltaPhiVarFn", ",", "\[IndentingNewLine]", 
    "DeltaPhiVar", ",", "\[IndentingNewLine]", "phiFn", ",", 
    "\[IndentingNewLine]", "dihedralHarmonicVarNames"}], "]"}], 
  ";"}]}], "Input",
 CellChangeTimes->{{3.952367425983066*^9, 3.952367555421852*^9}, {
  3.95236758915801*^9, 3.9523676657386723`*^9}, {3.952368034157886*^9, 
  3.9523680718003263`*^9}, {3.952368186579814*^9, 3.952368187104988*^9}, {
  3.952369295716834*^9, 3.95236930267612*^9}, {3.9523694272246*^9, 
  3.952369427602433*^9}},
 CellLabel->"In[43]:=",ExpressionUUID->"7ce95b23-e82c-4c40-9faa-b79443c43633"],

Cell[BoxData[
 RowBox[{"Length", "[", "dihedralHarmonicForceHessianRules", "]"}]], "Input",
 CellChangeTimes->{{3.9523755152998333`*^9, 3.9523755277711678`*^9}},
 CellLabel->"In[45]:=",ExpressionUUID->"3e8b5650-149b-4442-9e68-a35495317ad4"],

Cell[BoxData[
 RowBox[{
  RowBox[{"dihedralHarmonicRules", " ", "=", " ", 
   RowBox[{"Flatten", "[", " ", 
    RowBox[{"{", "\[IndentingNewLine]", 
     RowBox[{
     "dihedralHarmonicSetupRules", ",", "\[IndentingNewLine]", 
      "dihedralHarmonicEnergyRules", ",", "\[IndentingNewLine]", 
      "dihedralHarmonicForceHessianRules", ",", "\[IndentingNewLine]", 
      RowBox[{"CCode", "[", "\"\<SKIP_term: (void)0;\>\"", "]"}]}], "}"}], 
    "]"}]}], ";"}]], "Input",
 CellChangeTimes->{{3.952368160582329*^9, 3.9523681666217127`*^9}, {
  3.952368845636664*^9, 3.952368897280176*^9}, {3.9527492982489157`*^9, 
  3.952749299553296*^9}},
 CellLabel->"In[46]:=",ExpressionUUID->"ebc7700f-98aa-448d-bc3f-a53d4a2d9e04"],

Cell[BoxData[
 RowBox[{
  RowBox[{"dihedralHarmonicPack0", " ", "=", " ", 
   RowBox[{"{", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"Name", "\[Rule]", "\"\<DihedralHarmonicSpeedy\>\""}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"AdditionalCDeclares", "\[Rule]", "\"\<\>\""}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"EnergyFunction", "\[Rule]", "dihedralHarmonicEnergyFn"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"DerivativeVariables", "\[Rule]", "names"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{
     "HessianStructure", " ", "\[Rule]", " ", "dihedralHarmonicHessian"}], 
     ",", "\[IndentingNewLine]", 
     RowBox[{"Rules", "\[Rule]", "dihedralHarmonicRules"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"Input", "\[Rule]", "dihedralHarmonicInputs"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"Output", "\[Rule]", "dihedralHarmonicOutputs"}], ",", 
     "\[IndentingNewLine]", 
     RowBox[{"(*", 
      RowBox[{
       RowBox[{"ShortOptimizeCycles", "->", "True"}], ","}], "*)"}], 
     "\[IndentingNewLine]", 
     RowBox[{"(*", 
      RowBox[{
       RowBox[{"SkipTimesOptimize", "->", "True"}], ","}], "*)"}], 
     "\[IndentingNewLine]", 
     RowBox[{"SkipPlusOptimize", "->", "True"}]}], "\[IndentingNewLine]", 
    "}"}]}], ";"}]], "Input",
 CellChangeTimes->{{3.402880698908718*^9, 3.402880718144539*^9}, {
   3.405596810467842*^9, 3.405596826860199*^9}, {3.4056361468455667`*^9, 
   3.405636199543017*^9}, {3.407679917652017*^9, 3.407679931016655*^9}, 
   3.407763159631695*^9, {3.409333420940812*^9, 3.40933342773901*^9}, {
   3.410177864101145*^9, 3.410177865400531*^9}, {3.411066512729656*^9, 
   3.411066516761778*^9}, {3.41106697034696*^9, 3.411066970912582*^9}, {
   3.411073962997476*^9, 3.411073966619384*^9}, {3.411075733472687*^9, 
   3.411075737805724*^9}, {3.411076079199256*^9, 3.411076112584221*^9}, 
   3.411079373612965*^9, {3.9523689496208487`*^9, 3.952369017017926*^9}, {
   3.952592681842944*^9, 3.95259268948978*^9}, {3.952842818367293*^9, 
   3.952842825031496*^9}, {3.9528438720714607`*^9, 3.9528438828444242`*^9}, {
   3.95284476935154*^9, 3.9528447715548573`*^9}, {3.952849763281344*^9, 
   3.9528497741567173`*^9}, {3.9529353789979677`*^9, 3.952935382828443*^9}, {
   3.9529567652322893`*^9, 3.95295676612022*^9}, {3.9529689032600946`*^9, 
   3.952968906120604*^9}, {3.9529770072504807`*^9, 3.952977015128125*^9}, {
   3.952985635099257*^9, 3.952985638159172*^9}},
 CellLabel->"In[47]:=",ExpressionUUID->"5c35519d-f10a-4464-a87b-947091ce09ee"],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{
   RowBox[{"dihedralHarmonicPack", " ", "=", " ", 
    RowBox[{"packOptimize", "[", "dihedralHarmonicPack0", "]"}]}], ";"}], 
  "*)"}]], "Input",
 CellChangeTimes->{{3.952369067607729*^9, 3.9523690763868837`*^9}, {
  3.952376086877322*^9, 3.952376089938352*^9}, {3.952382801853137*^9, 
  3.9523828036970863`*^9}, {3.952940878281062*^9, 3.9529408820364933`*^9}},
 CellLabel->"In[48]:=",ExpressionUUID->"101db157-22e0-4823-8a0b-b1fcc3602649"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"dihedralHarmonicPack", " ", "=", " ", 
   RowBox[{"packOptimize", "[", "dihedralHarmonicPack0", "]"}]}], 
  ";"}]], "Input",
 CellChangeTimes->{{3.9529472026068993`*^9, 3.952947219137217*^9}},
 CellLabel->"In[49]:=",ExpressionUUID->"540ebd5d-9c4c-4fdc-85b5-7fb4e260c9c7"],

Cell[BoxData[
 RowBox[{
  RowBox[{"rules", " ", "=", " ", 
   RowBox[{"Rules", "/.", "dihedralHarmonicPack"}]}], ";"}]], "Input",
 CellChangeTimes->{{3.952966188649509*^9, 3.952966214595786*^9}},
 CellLabel->"In[50]:=",ExpressionUUID->"3c71f69a-3850-490f-a199-b8ddd289d3e6"],

Cell[BoxData[
 RowBox[{"Length", "[", "rules", "]"}]], "Input",
 CellChangeTimes->{{3.95296689471618*^9, 3.952966896993132*^9}},
 CellLabel->"In[51]:=",ExpressionUUID->"40e2f135-57f0-45a7-a91a-458f73892365"],

Cell[BoxData[
 RowBox[{"EchoTiming", "[", 
  RowBox[{"associativeMostCommonPairMap", "[", 
   RowBox[{"rules", ",", "Times"}], "]"}], "]"}]], "Input",
 CellChangeTimes->{{3.952966224938822*^9, 3.952966236348068*^9}, {
   3.952966817651689*^9, 3.952966888320166*^9}, 3.952967124358678*^9, {
   3.9529672240939503`*^9, 3.9529672381925507`*^9}},
 CellLabel->"In[52]:=",ExpressionUUID->"5b29f4f1-176e-44bc-85b2-9994907a1ac5"],

Cell[BoxData[
 RowBox[{"rules", " ", "=", " ", 
  RowBox[{"Rules", "/.", "dihedralHarmonicPack"}]}]], "Input",
 CellChangeTimes->{{3.95298291863081*^9, 3.952982930464497*^9}},
 CellLabel->"In[53]:=",ExpressionUUID->"dc27fb53-b9e4-410f-9ff8-97b83ccb2222"],

Cell[BoxData[
 RowBox[{
  RowBox[{"prods", " ", "=", " ", 
   RowBox[{"associativeExtractComplexTerms", "[", 
    RowBox[{"rules", ",", "Times"}], "]"}]}], ";"}]], "Input",
 CellChangeTimes->{{3.952982933135077*^9, 3.952982969291645*^9}},
 CellLabel->"In[54]:=",ExpressionUUID->"28ec59f3-ad9b-46ba-a8e0-05e4b875e151"],

Cell[BoxData[
 RowBox[{"Length", "[", "prods", "]"}]], "Input",
 CellChangeTimes->{{3.952982973106525*^9, 3.952982978303812*^9}},
 CellLabel->"In[55]:=",ExpressionUUID->"b64e72e6-0368-4a59-b0c8-18bc428bd260"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"index", "=", 
    RowBox[{"EchoTiming", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"Reap", "[", 
        RowBox[{
         RowBox[{"Do", "[", 
          RowBox[{
           RowBox[{
            RowBox[{"symbols", "=", 
             RowBox[{"List", "@@", 
              RowBox[{"rule", "[", 
               RowBox[{"[", "1", "]"}], "]"}]}]}], ";", "\[IndentingNewLine]", 
            RowBox[{"Scan", "[", 
             RowBox[{
              RowBox[{
               RowBox[{"Sow", "[", 
                RowBox[{"rule", ",", "#"}], "]"}], "&"}], ",", "symbols"}], 
             "]"}]}], ",", 
           RowBox[{"{", 
            RowBox[{"rule", ",", "prods"}], "}"}]}], "]"}], ",", "_", ",", 
         "Rule"}], "]"}], "[", 
       RowBox[{"[", "2", "]"}], "]"}], ",", "\"\<ReapSow index\>\""}], 
     "]"}]}], ";"}], "\n"}]], "Input",
 CellChangeTimes->{{3.95298299349371*^9, 3.9529829935025663`*^9}, 
   3.952983044846304*^9, {3.952983085926272*^9, 3.952983086681011*^9}},
 CellLabel->"In[56]:=",ExpressionUUID->"2638ef08-cc88-456b-95b2-3db19df43c1f"],

Cell[BoxData[
 RowBox[{"ass", " ", "=", " ", 
  RowBox[{"Association", "[", "index", "]"}]}]], "Input",
 CellChangeTimes->{{3.9529829995743837`*^9, 3.9529830002330103`*^9}, {
  3.9529830995417013`*^9, 3.9529831289305897`*^9}},
 CellLabel->"In[57]:=",ExpressionUUID->"c64b00ac-fabd-44cc-8abd-6c81987587e1"],

Cell[BoxData[
 RowBox[{"ass", "[", 
  RowBox[{"[", 
   RowBox[{"Key", "[", "tx15", "]"}], "]"}], "]"}]], "Input",
 CellChangeTimes->{{3.9529831792393913`*^9, 3.952983184093041*^9}},
 CellLabel->"In[58]:=",ExpressionUUID->"3d4ed4b0-3e9c-4b4f-a14d-cbfed8a420cf"],

Cell[BoxData[
 RowBox[{"outputs", " ", "=", " ", 
  RowBox[{"{", "}"}]}]], "Input",
 CellChangeTimes->{{3.9529774949681664`*^9, 3.952977496480213*^9}},
 CellLabel->"In[59]:=",ExpressionUUID->"783b98c1-bc97-4457-97b1-0e16dc2ff853"],

Cell[BoxData[
 RowBox[{"associativesSimplify", "[", 
  RowBox[{"dihedralHarmonicPack", ",", "Times", ",", "outputs"}], 
  "]"}]], "Input",
 CellChangeTimes->{{3.9529774660673122`*^9, 3.9529775105459003`*^9}},
 CellLabel->"In[60]:=",ExpressionUUID->"e4515d97-66e4-4b22-b7ef-0a27b7731cb1"]
},
Evaluator->"Foo",
CellGrouping->Manual,
WindowSize->{1533, 1067},
WindowMargins->{{-2696, Automatic}, {Automatic, 88}},
FrontEndVersion->"13.3 for Mac OS X ARM (64-bit) (July 24, 2023)",
StyleDefinitions->"Default.nb",
ExpressionUUID->"953b3b28-3dbc-48a7-818e-7bd57abcc6b4"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[580, 22, 293, 6, 98, "Title",ExpressionUUID->"59ad4f77-f93b-4129-8c5d-5c71330c705c",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[876, 30, 628, 11, 30, "Input",ExpressionUUID->"2302b7fc-c3e0-473a-9944-bd8907480542",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[1507, 43, 352, 6, 52, "Code",ExpressionUUID->"571b4188-8844-47fd-a5e9-9e8bc63da377",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[1862, 51, 841, 13, 52, "Code",ExpressionUUID->"13260320-554d-4f6a-97fe-a104adcbae00",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[2706, 66, 425, 7, 53, "Subtitle",ExpressionUUID->"f3c71499-dcd8-4a60-8091-9938479e0665",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[3134, 75, 289, 5, 52, "Code",ExpressionUUID->"29cf2b5b-9ced-479c-85ae-c0a9be6c68e4",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[3426, 82, 381, 8, 52, "Code",ExpressionUUID->"8c98e9fb-8d7f-4365-9e46-d53e1da154f8",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[3810, 92, 464, 7, 52, "Code",ExpressionUUID->"420f0456-2a09-4d78-bd8d-6969828b5b21",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}]
}, Open  ]],
Cell[CellGroupData[{
Cell[4311, 104, 431, 10, 69, "Chapter",ExpressionUUID->"1f1aa4cc-2a9c-461d-93ec-97c2242c265e",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[4745, 116, 757, 13, 30, "Input",ExpressionUUID->"eda73ac1-ca89-4acb-93f8-a2bcaef0ea63",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[5505, 131, 1028, 26, 94, "Input",ExpressionUUID->"ab8ce859-bbb6-4364-b2a0-760d04215a00",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[6536, 159, 350, 7, 30, "Input",ExpressionUUID->"f4046a73-8241-4cc3-92b5-b0168510db9b",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[6889, 168, 349, 7, 30, "Input",ExpressionUUID->"eda20186-450b-48df-880f-677c611f8848",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[7241, 177, 2821, 67, 346, "Input",ExpressionUUID->"9d21da8e-a00b-4bbf-819a-d0b689cd5914",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[10065, 246, 447, 10, 30, "Input",ExpressionUUID->"84b668a5-368f-441a-98ba-2cff9272ae57",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[10515, 258, 595, 12, 30, "Input",ExpressionUUID->"115dfce1-7a62-4a15-8b38-98017ab13dfc",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[11113, 272, 799, 16, 30, "Input",ExpressionUUID->"f7032a5f-8737-4c9b-9eca-017969e08d46",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[11915, 290, 1251, 21, 46, "Input",ExpressionUUID->"c3d70802-8e89-4cfd-aeb7-755b7ecd4da5",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[13169, 313, 454, 8, 30, "Input",ExpressionUUID->"d4b06aff-2869-4380-8cb9-fe70fe3b34b8",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[13626, 323, 1122, 25, 136, "Input",ExpressionUUID->"8b9e4bcf-6850-4d82-bade-26d978769d42",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[14751, 350, 844, 22, 30, "Input",ExpressionUUID->"34e7ce02-3852-496d-9c58-948a1952826d",
 CellGroupingRules->{"GroupTogetherGrouping", 10000.}],
Cell[15598, 374, 328, 7, 30, "Input",ExpressionUUID->"266cf776-52f1-4a43-b178-bdf8459ae224"],
Cell[15929, 383, 417, 9, 30, "Input",ExpressionUUID->"590d8c9e-4792-413c-9810-3bb0c7f7e9a1"],
Cell[16349, 394, 158, 3, 54, "Subsection",ExpressionUUID->"6fb2043b-fed8-4159-9046-f2a26ced5774"],
Cell[16510, 399, 2031, 46, 300, "Input",ExpressionUUID->"42d21f14-444c-49d2-9590-1099b09c7d89"],
Cell[18544, 447, 2134, 52, 178, "Input",ExpressionUUID->"a18728c5-d7d6-43e3-8ba0-e33e1aeffaa8"],
Cell[20681, 501, 1761, 38, 136, "Input",ExpressionUUID->"db4a22fd-feda-40c4-942d-77cc9d3ce06a"],
Cell[22445, 541, 328, 5, 30, "Input",ExpressionUUID->"9888ee28-cae6-446b-a9ed-7bb0a5b91811"],
Cell[22776, 548, 293, 6, 30, "Input",ExpressionUUID->"c1775d99-a2a5-4151-95e4-b6accbfb9c97"],
Cell[23072, 556, 468, 10, 30, "Input",ExpressionUUID->"d45175be-4779-4d6a-ba17-88e005b2e201"],
Cell[23543, 568, 2925, 54, 136, "Input",ExpressionUUID->"fa831a18-e7d4-42d3-86d1-85c6f7869e16"],
Cell[26471, 624, 1103, 20, 220, "Input",ExpressionUUID->"7ce95b23-e82c-4c40-9faa-b79443c43633"],
Cell[27577, 646, 240, 3, 30, "Input",ExpressionUUID->"3e8b5650-149b-4442-9e68-a35495317ad4"],
Cell[27820, 651, 717, 14, 115, "Input",ExpressionUUID->"ebc7700f-98aa-448d-bc3f-a53d4a2d9e04"],
Cell[28540, 667, 2540, 48, 283, "Input",ExpressionUUID->"5c35519d-f10a-4464-a87b-947091ce09ee"],
Cell[31083, 717, 481, 9, 30, "Input",ExpressionUUID->"101db157-22e0-4823-8a0b-b1fcc3602649"]
}, Open  ]],
Cell[31579, 729, 306, 6, 30, "Input",ExpressionUUID->"540ebd5d-9c4c-4fdc-85b5-7fb4e260c9c7"],
Cell[31888, 737, 274, 5, 30, "Input",ExpressionUUID->"3c71f69a-3850-490f-a199-b8ddd289d3e6"],
Cell[32165, 744, 207, 3, 30, "Input",ExpressionUUID->"40e2f135-57f0-45a7-a91a-458f73892365"],
Cell[32375, 749, 421, 7, 30, "Input",ExpressionUUID->"5b29f4f1-176e-44bc-85b2-9994907a1ac5"],
Cell[32799, 758, 254, 4, 30, "Input",ExpressionUUID->"dc27fb53-b9e4-410f-9ff8-97b83ccb2222"],
Cell[33056, 764, 317, 6, 30, "Input",ExpressionUUID->"28ec59f3-ad9b-46ba-a8e0-05e4b875e151"],
Cell[33376, 772, 208, 3, 30, "Input",ExpressionUUID->"b64e72e6-0368-4a59-b0c8-18bc428bd260"],
Cell[33587, 777, 1106, 29, 73, "Input",ExpressionUUID->"2638ef08-cc88-456b-95b2-3db19df43c1f"],
Cell[34696, 808, 305, 5, 30, "Input",ExpressionUUID->"c64b00ac-fabd-44cc-8abd-6c81987587e1"],
Cell[35004, 815, 260, 5, 30, "Input",ExpressionUUID->"3d4ed4b0-3e9c-4b4f-a14d-cbfed8a420cf"],
Cell[35267, 822, 230, 4, 30, "Input",ExpressionUUID->"783b98c1-bc97-4457-97b1-0e16dc2ff853"],
Cell[35500, 828, 287, 5, 30, "Input",ExpressionUUID->"e4515d97-66e4-4b22-b7ef-0a27b7731cb1"]
}
]
*)

