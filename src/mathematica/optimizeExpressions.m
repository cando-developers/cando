(* ::Package:: *)

(* ::Subtitle:: *)
(*Generate optimized energy, gradient, and Hessian expressions for molecular mechanics.*)


(* ::Subsubtitle:: *)
(*collectTermsForList( equations, terms ) -> terms*)
(**)
(*Collects terms for a list of equations and returns a list of terms that calculates the equations in a more efficient manner.*)
(*All sub-expressions are evaluated once, put into temporary variables and used over and over again to calculate all of the subsequent equations.*)


(* ::Subsubtitle:: *)
(*Calculate the gradient and Hessian using a list of variables.*)
(*Only the upper triangular elements of the Hessian are calculated.*)


(* ::Text:: *)
(*Turn off those annoying spelling warnings*)


BeginPackage["optimizeExpressions`"]


Off[General::spell];
Off[General::spell1];


oneForceName[str_] := "f"<>ToString[str];
oneDiagonalHessianName[str_] := "dh"<>ToString[str]<>ToString[str];
oneOffDiagonalHessianName[str_] := "oh"<>ToString[str];


firstDeriv[fnE_,c_,fnC_,v_] := D[fnE,c]D[fnC,v]


secondDeriv[fnE_,c_,fnC_,v1_,v2_] := D[D[fnE,c],c]D[fnC,v2]D[fnC,v1]+D[fnE,c]D[D[fnC,v2],v1]


buildGradient[prefix_,fnE_,c_,fnC_,names_] := Block[{deriv,i},
	Table[firstDeriv[fnE,c,fnC,names[[i]]]->Symbol[ToString[prefix]<>ToString[names[[i]]]],{i,Length[names]}]
];


buildHessian[prefix_,fnE_,c_,fnC_,names_] := Block[{hessian,i,j},
	hessian = {};
	For[i=1,i<=Length[names], i++,
		For[j=i,j<=Length[names], j++,
			hessian = Append[hessian,secondDeriv[fnE,c,fnC,names[[i]],names[[j]]]->Symbol[ToString[prefix]<>ToString[names[[i]]]<>ToString[names[[j]]]]];
		];
	];
	hessian
];


(* ::Subsection:: *)
(*selectRules selects out only rules from the list of expressions.*)


selectRules[exps_] := Select[exps,Head[#]===Rule&];


packSelectRules[pack_] := Select[Rules/.pack,Head[#]==Rule&]


simplifyEquation[eq_,terms_]:= Module[{ne},
ne = eq/.selectRules[terms];
If[eq==ne,Return[ne];];
ne = simplifyEquation[ne,terms];
Return[ne];
];


collectTermsForOneLevel[fAndTerms_] := Module[{newExpression,lev,ts,nt,expr, nf,i},
	nf = fAndTerms[[1]];
	nt = fAndTerms[[2]];
	(*nf = nf/.nt;*)
	lev = Depth[nf];
	ts = Union[Level[nf,{lev-2}]](*/.nt*);
	For[i=1,i<=Length[ts],i++,
		expr = ts[[i]];
		If[!((Head[expr]=== Integer) ||
			(Head[expr]=== Real) ||
			(Head[expr]=== Symbol) ||
			(Head[expr] === Rational)),
			(* then *) map[expr] = Symbol["tx"<>ToString[uniqueVariableIndex]];(*Variable definition*)
						uniqueVariableIndex = uniqueVariableIndex+1;
						newExpression = expr->map[expr];
						nt = Append[nt,newExpression];
						PrintTemporary["Newish expression: "<>ToString[newExpression]]; (*Comment out*)
						nf =  nf/.(expr -> map[expr]);
		];
	];
	Return[{nf,nt}];
];


collectTermsForOneExpression[name_,fn_, terms_]:=Module[{nf,nt,res},
	(* PrintTemporary["name=",name," fAndTerms=",fAndTerms]; *)
	nt=terms;
	nf=simplifyEquation[fn,nt];
	While[Depth[nf]>2,
		res=collectTermsForOneLevel[{nf,nt}];
		nf=res[[1]];
		nt=res[[2]];
	];
	nt=Append[nt,nf->name];
	Return[nt];
];


(* ::Subsection:: *)
(*Remove all of the trivial rules from the list of rules.*)
(*Rules like (a->b) replace any occurance of b with a.*)


trivialRules[rules_,outs_] := Block[{triv},
	triv = Select[rules,Head[#]===Rule && Head[#[[1]]]===Symbol && Count[outs,#[[1]]]==0&];
	Print["trivialRules>> triv = "<>ToString[triv]];
	Print["trivialRules>> outs = "<>ToString[outs]];
	triv
];


oldNonTrivialRules[rules_,trivial_] := Block[{nontriv},
	Select[rules,Count[trivial,#]==0&];
	nontriv
];


nonTrivialRules[rules_,outs_] := Block[{nontriv},
	nontriv = Select[rules,Not[Head[#]===Rule && Head[#[[1]]]===Symbol && Count[outs,#[[1]]]==0]&];
	nontriv
];


eliminateTrivialRules[packet_,outputs_] := Module[{rules,trivRules,nonTrivRules,leftRules},
	rules = Rules/.packet;
	trivRules = trivialRules[rules,outputs];
	If[Length[trivRules]>0,
		(* then *) nonTrivRules = rules; (* nonTrivialRules[rules,trivRules]; *)
				(*nonTrivRules = nonTrivRules//.trivRules;*)
				Print["Trivial rules are being removed: "<>ToString[trivRules//MatrixForm]];
				Print["After removal: " <> ToString[nonTrivRules//MatrixForm]];
				Return[{Rules->nonTrivRules,Changes->True}];
		(* else *),
				Print["There were no trivial rules"];
	];
	{Rules->rules,Changes->False}
];


(* ::Subsection:: *)
(*For a list of equations assigned to variables as rules, collect all sub-expressions that are common and assign them to temporary variables.*)
(*      Return the new list of equations assigned to variables.*)
(*      Use a packet which is a list of rules to bring the data into and out of this function.*)
(*      The packet should have a rule: "Rules->{xxxx}" coming in*)
(*      and will return with a rule: "Rules->" and "Changes->True/False" if there were changes to the rules*)


collectTerms[packet_,outputs_] := Module[{changes,pack,equations,newnt,result, nt,i,numberOfNewTerms},
equations = Rules/.packet;
nt = {};
changes = False;
For[i=1,i<=Length[equations],i++,
If[ Head[equations[[i]]]===Rule, (*If this is a rule then collect the sub-expressions, if not then just add it to the list of equations *)
(*PrintTemporary["Collecting terms for rule: ",i, " ", equations[[i]][[2]]];*)
newnt = collectTermsForOneExpression[equations[[i]][[2]],equations[[i]][[1]],nt];
numberOfNewTerms = Length[newnt]-Length[nt]-1; (* One term is always added from the equations list *)
If[numberOfNewTerms!=0,
changes = True;
PrintTemporary["   collected "<>ToString[numberOfNewTerms]<>" terms"]
];
nt = newnt;
,
(*Print["Appending non-term"<>equations[[i]]];*)
nt = Append[nt,equations[[i]]];
];
(* Print[nt]; *) (*Comment out*)
];
pack = { Rules->nt};
Print["eliminateTrivialRules"];
result =eliminateTrivialRules[pack,outputs];
changes = changes||(Changes/.result);
result ={Rules -> (Rules/.result),Changes->changes};
If[changes,PrintTemporary["CHANGES MADE, RERUN"]];
result
];


(* ::Section:: *)
(*simplify products and sums*)
(**)
(*Mathematica represents products with the Times[x,y,z,...] function.  These need to be simplified because they can contain pair products that are used in other expressions.  The following code identifies common pair products, creates a new term/variable for it, inserts it into the term list just before the first term that would require them and substitutes every term that requires it with the new variable.*)
(**)
(*The -simplifyProducts- function keeps operating on the terms list until all pair products that are evaluated more than once are simplified.*)


associativePrintComplexTerms[terms_,type_] := Module[{allVars,allTerms,i,j,ot,allNt,nt},
allTerms = Select[selectRules[terms],Head[#[[1]]]===type && Length[#[[1]]]>2 &];
allVars = {};
For[i=1,i<=Length[allTerms],i++,
ot=Table[allTerms[[i]][[1]][[j]],{j,1,Length[allTerms[[i]][[1]]]}];
allVars = Union[allVars,ot];
];
If[Length[allVars]==0,Return["All Done"]];
allNt={};
For[i=1,i<=Length[allTerms],i++,
nt = {};
For[j=1,j<=Length[allVars],j++,
If[MemberQ[allTerms[[i]][[1]],allVars[[j]]],nt = Append[nt,allVars[[j]]];,nt=Append[nt," "];];
];
allNt = Append[allNt,nt];
];
Return[Sort[allNt]];
];


associativeTermContains[assoc_,vars_] := Module[{tl,i,res},
(*Return True if the associative object contains all of the vars.*)
tl = Table[assoc[[i]],{i,1,Length[assoc]}]; (*Accumulate all of the variables in the associative term*)
If[Union[tl,vars]==tl,res=True,res=False]; (*Return True if vars is a subset of t1*)
Return[res]
];


associativeSubstitute[rule_,vars_,newVar_,type_] := Module[{tnew,newTimes,leftSide,rightSide,tl,i,newTerm},
(*Print["productSubstitute term=",term," vars=",vars," newVar=",newVar];*)
leftSide = rule[[1]];
rightSide = rule[[2]];
tl = Level[leftSide,1];
tnew = Append[Complement[tl,vars],newVar];
newTimes = Apply[type,tnew];
(*Print["associativeSubstitute newTimes=",newTimes];*)
newTerm = newTimes->rightSide;
(*Print["associativeSubstitute newTerm=", newTerm];*)
Return[newTerm];
];


associativeFirstUse[rules_,vars_,type_] := Module[{i,use},
(*Print["associativeFirstUse:: rules=",rules];*)
(*Print["associativeFirstUse:: vars=",vars];*)
(*Print["associativeFirstUse:: type=",type];*)
For[i=1,i<=Length[rules],i++,
If[Head[rules[[i]]]===Rule,
If[Head[rules[[i]][[1]]]===type,
If[associativeTermContains[rules[[i]][[1]],vars],Return[i]];
];
];
];
Return[-1]];


associativeSimplifyPair[rules_,vars_,type_] := Module[{pos,fu,i,newRules,subs, oneRule,insTerm,insertedTerms, newTimes, newVar},
(*Create a new variable for the pair product and add its definition to the term list and then substitute all matching pair products with the new definition*)
(*Print["associativeSimplifyPair:: starting"];*)
pos = associativeFirstUse[rules,vars,type];
(*Print["First use=",pos];*)
newVar = Symbol["tzz"<>ToString[uniqueVariableIndex]];
uniqueVariableIndex = uniqueVariableIndex + 1;
(*Print["New var= ",newVar ];*)
If[pos<0,Throw["terms not found: "<>ToString[vars]]];
newRules={};
subs = 0;
For[i=1,i<=Length[rules],i++,
oneRule = rules[[i]];
If[Head[oneRule]===Rule,
If[Head[oneRule[[1]]]===type,
If[associativeTermContains[oneRule[[1]],vars],
(*Print["associativeSimplifyPair:: oneRule=", oneRule];*)
oneRule=associativeSubstitute[oneRule,vars,newVar,type];
(*Print["associativeSimplifyPair new oneRule=", oneRule];*)
subs++;
];
];
];
AppendTo[newRules,oneRule];
];
insTerm = Apply[type,vars]->newVar;
insertedTerms = Insert[newRules,insTerm,pos];
PrintTemporary["Substitute term: "<>ToString[CForm[insTerm[[2]]]] <>"="
				<>ToString[CForm[insTerm[[1]]]]<>" at position:"
				<>ToString[pos]<>" substitutions: "<>ToString[subs]];
Return[insertedTerms];
];


associativeExtractComplexTerms[terms_,type_] := Select[selectRules[terms],Head[#[[1]]]===type&&Length[#[[1]]]>2 &];


associativeExtractFactors[rules_] := Table[Table[rules[[i]][[1]][[j]],{j,1,Length[rules[[i]][[1]]]}],{i,1,Length[rules]}];


associativePossiblePairs[listOfFactors_] := Module[{a,i,j,factors,pp},
pp = {};
For[a=1,a<=Length[listOfFactors],a++,
factors = Sort[listOfFactors[[a]]];
For[i=1,i<=Length[factors]-1,i++,
For[j=i+1,j<=Length[factors],j++,
AppendTo[pp,{factors[[i]],factors[[j]]}];
];
];
];
Return[Union[pp]];
];


associativeMostCommonPair[terms_,type_] := Module[{prods, all,pair,pairs,seen,trialPairs,i,j},
prods = associativeExtractComplexTerms[terms,type];
If[Length[prods]==0,Return[{}]];
(*PrintTemporary["Prods=",prods];*)
all = associativeExtractFactors[prods];
trialPairs = associativePossiblePairs[all];
(*PrintTemporary["Possible pairs=",trialPairs ];*)
pairs = {};
For[i=1,i<Length[trialPairs],i++,
pair = trialPairs[[i]];
(*PrintTemporary["Pair = ", pair];*)
seen = Cases[prods,x_/;MemberQ[x[[1]],pair[[1]]]&&MemberQ[x[[1]],pair[[2]]]];
(*PrintTemporary[ "seen=", seen];*)
pairs = Append[pairs,{Length[seen],pair,seen}];
];
Return[Sort[pairs][[-1]]];
];


associativesSimplify[packet_,type_,outputs_]:=Module[{terms,newTerms,done,common},
(*Right now simplifyProducts doesn't set Changes properly when a pair product is simplified.*)
	terms = Rules/.packet;
	newTerms = terms;
	done = False;
	While[ !done,
		common = associativeMostCommonPair[newTerms,type];
		If[Length[common]==0,done= True,
			(*PrintTemporary["Most common: ",common[[2]], " number=",common[[1]]];*)
			If[common[[1]]>1,
				newTerms=associativeSimplifyPair[newTerms,common[[2]],type];
			,
				done=True;
			];
			(*PrintTemporary[printComplexProductTerms[newTerms]//MatrixForm];*)
		];
	];
	PrintTemporary["associativeSimplify>>Got to eliminateTrivialRules"];
	newTerms = eliminateTrivialRules[{Rules->newTerms},outputs];
	(*PrintTemporary[newTerms];*)
	PrintTemporary["associativeSimplify>>Came out of eliminateTrivialRules"];
	{Rules->(Rules/.newTerms),Changes->(Changes/.newTerms)}
];


timesSimplify[packet_,outs_] := associativesSimplify[packet,Times,outs];
plusSimplify[packet_,outs_] := associativesSimplify[packet,Plus,outs];


(* ::Section:: *)
(*Write out the list of terms in a form that can be cut and pasted into C code.*)


variableList[terms_,outputs_] := Block[
	{rules,variables,on,res},
	rules = selectRules[terms];
	on = ToString/@outputs;
	res = Union[Table[ToString[rules[[i]][[2]]],{i,Length[rules]}],on];
	Return[res];
];


cFormVariableDeclarations[terms_,outputs_,additionalCDeclares_] := Module[{u,i},
	u = variableList[terms,outputs];
	decls = "";
	For[i=1,i<=Length[u],i++,decls = decls<>"DECLARE_FLOAT("<>ToString[u[[i]]]<>");\n"];
	decls = decls<>additionalCDeclares<>"\n";
	Return[decls];
];


toC[Rule[lhs_,rhs_],idx_] := Block[{}, "\t "<>ToString[rhs//CForm] <> " = " <> ToString[lhs//CForm]<>"; \t\t/* rule "<>ToString[idx]<>" */\n"]


toC[CCode[str_String],idx_] := Block[{},"\t"<>str<>"\n"];


toC[MathCode[exp_],idx_] := Block[{},"/*"<>ToString[exp]<>"*/\n"];


cFormCode[terms_] := Block[{i,str},
	str = "// C-code\n";
	For[i=1,i<=Length[terms],i++,
		str = str <> toC[terms[[i]],i];
	];
	Return [str];
];


cFormListOfTerms[terms_,outputs_,additionalCDeclares_] := Module[{i,str},
	str = "#ifdef TERM_DECLARES //[\n";
	str = str<>cFormVariableDeclarations[terms,outputs,additionalCDeclares]<>"\n";
	str = str<>"#endif //TERM_DECLARES ]\n";
	str = str<>"#ifdef TERM_CODE //[\n";
	For[i=1,i<=Length[terms],i++,
		str = str <> toC[terms[[i]],i];
		(* If[Head[terms[[i]]]===Rule,
			(*then*) str = str <> "\t"<>ToString[terms[[i]][[2]]//CForm] <> " = " <> ToString[terms[[i]][[1]]//CForm]<>"; \t\t/* rule "<>ToString[i]<>" */\n";,
			(*else*) str = str <> "\t"<>ToString[terms[[i]]]<>"\n";
		]; *)
	];
	str = str<>"#endif //TERM_CODE ]\n";
	Return[str];
];


cFormWriteToFiles[rules_,pack_]:=Block[
	{fileName,st,declares,code},
	fileName = "_"<>(Name/.pack)<>"_termDeclares.cc";
Print["Writing declares to file: ", fileName ];
	str = "// declare variables\n";
    str = str<>cFormVariableDeclarations[rules,(Output/.pack),(AdditionalCDeclares/.pack)]<>"\n";
	st = OpenWrite[fileName];
	WriteString[st,"#if USE_EXPLICIT_DECLARES\n"];
	WriteString[st,str];
	WriteString[st,"#endif\n"];
	Close[st];
	fileName = "_"<>(Name/.pack)<>"_termCode.cc";
	code = cFormCode[rules];
	Print["Writing code to file: ", fileName];
	st = OpenWrite[fileName];
	WriteString[st,code];
	Close[st];
];


(* ::Section:: *)
(*Generate Mathematica block that evaluates the pack*)


mathComment[idx_] := " (*Rule "<>ToString[idx]<>"*)\n";


toMath[Rule[lhs_,rhs_],idx_] := "\t"<>ToString[rhs] <>"="<>ToString[lhs]<>";"<>mathComment[idx];


toMath[Rule[power2[val_],rhs_],idx_] := "\t"<>ToString[rhs]<>"= mathPower2["<>ToString[val]<>"];"<>mathComment[idx];


toMath[Rule[mysqrt[val_],rhs_],idx_] := "\t"<>ToString[rhs]<>"= mathSqrt["<>ToString[val]<>"];"<>mathComment[idx];


toMath[Rule[reciprocal[val_],rhs_],idx_] := "\t"<>ToString[rhs]<>"= mathReciprocal["<>ToString[val]<>"];"<>mathComment[idx];


toMath[CCode[str_],idx_] := "\t"<>"(*"<>ToString[str]<>"*)\n";


toMath[MathCode[stuff_],idx_] := "\t"<>ToString[stuff]<>"\n";


packAsMathematicaBlock[rules_,outputs_] := Block[{str,vars},
	str = "Block[{";
	vars = variableList[rules,outputs];
	str = str <> StringJoin[Table[vars[[i]]<>",",{i,Length[vars]}]];
	str = str <> "xxxDummy},\n";
	str = StringJoin[{str,Table[toMath[rules[[i]],i],{i,Length[rules]}]}];
	str = str<>ToString[outputs];
	str = str<>"]\n";
	str
];


(* ::Section:: *)
(*Generate Mathematica code to test a list of terms to see if it generates the right energy/gradient/Hessian*)


(* ::Text:: *)
(*AppendTo[br, mathematicaSinNPhi[IN, sinPhi, cosPhi] -> SinNPhi];*)
(*AppendTo[br, mathematicaCosNPhi[IN, sinPhi, cosPhi] -> CosNPhi];*)


testCaseExpressions = "
Begin[\"`Private`\"];
oneHalf = 0.5;
two = 2.0;
twelve=12.0;
six=6.0;
EraseLinearDihedral=1.0;
power2[x_] = x*x;
mysqrt[x_] = Sqrt[x];
reciprocal[x_] = 1.0/x;
mathematicaSinNPhi[1,sinPhi_,cosPhi_] = sinPhi;
mathematicaCosNPhi[1,sinPhi_,cosPhi_] = cosPhi;
mathematicaSinNPhi[n_,sinPhi_,cosPhi_] = cosPhi*mathematicaSinNPhi[n-1,sinPhi,cosPhi]+sinPhi*mathematicaCosNPhi[n-1,sinPhi,cosPhi];
mathematicaCosNPhi[n_,sinPhi_,cosPhi_] = cosPhi*mathematicaCosNPhi[n-1,sinPhi,cosPhi]-sinPhi*mathematicaSinNPhi[n-1,sinPhi,cosPhi];
";
testCaseExpressionsClear = "Clear[oneHalf,two,twelve,six,power2,mysqrt,reciprocal,mathematicaSinNPhi,mathematicaCosNPhi];";


testCaseBlock[interaction_,pack_] := Module[{results,args,rules,body,vars,v,str,termRules,i,result},
rules = Rules/.pack;
args = Args/.interaction;
results=Results/.interaction;
vars={};
(*terminals = {};*)
body = testCaseExpressions;
For[i=1,i<=Length[args],i++,
body = body <>ToString[args[[i]][[1]]]<>"="<>ToString[args[[i]][[2]]//InputForm]<>";\n";
(*body = body<>"Print[\""<>ToString[trials[[i]]]<>"\"];\n";*)
AppendTo[vars,args[[i]][[1]]];
];
(*For[i=1,i<=Length[results],i++,
(*body = body<>"Print[\""<>ToString[trials[[i]]]<>"\"];\n";*)
AppendTo[vars,results[[i]][[1]]];
];*)
termRules = selectRules[rules];
For[i=1,i<=Length[termRules],i++,
v = termRules[[i]][[2]];
body = body<>ToString[v]<>"="<>ToString[termRules[[i]][[1]]//InputForm]<>";\n";
AppendTo[vars,v];
(*If[StringTake[ToString[v],1]=="t",AppendTo[terminals,v]];*)
];
(*For[i=1,i<=Length[terminals],i++,
body=body<>"Print[\""<>ToString[terminals[[i]]]<>"=\","<>ToString[terminals[[i]]]<>"];\n";
];*)
body = body <> "{";
For[i=1,i<=Length[vars],i++,
(*body=body<>ToString[vars[[i]]]<>"->"<>"\""<>ToString[vars[[i]]]<>"\"";*)
body=body<>"\""<>ToString[vars[[i]]]<>"\""<>"->"<>ToString[vars[[i]]];
If[i<Length[vars],body=body<>","];
];
body=body<>"}";
str ="Block["<>ToString[vars]<>",\n"<>body<>"]";
str
];


testCaseBlockEvaluate[tb_] := Module[{result},
	result = ToExpression[tb];
	ToExpression[testCaseExpressionsClear];
	Table[Symbol[result[[i]][[1]]]->result[[i]][[2]],{i,Length[result]}]
];


testCaseEvaluate[args_,pack_] := Module[{tb,result},
	tb = testCaseBlock[args,pack];
	testCaseBlockEvaluate[tb]
];


testCaseCompareBlock[interaction_,tb_]:=Module[
	{mathematicaResults,interResults,comparison,res,ires,mres,var,i},
	mathematicaResults = testCaseBlockEvaluate[tb];
	interResults = Results/.interaction;
	PrintTemporary["interResults="<>ToString[interResults]];
	comparison = {{"var","amber/cando","mathematica","difference","%difference"}};
	For[i=1,i<=Length[interResults],i++,
		var = interResults[[i]][[1]];
		ires = interResults[[i]][[2]];
		mres = var/.mathematicaResults;
		AppendTo[comparison,{var,ires,mres,ires-mres,100*Abs[ires-mres]/(Abs[ires+mres]/2.0)}];
	];
	comparison
];


testCaseCreateCompiledFunction[interaction_,pack_] := Block[{},
	ToExpression[MathCode/.pack]
];


testCaseCompareCompiledFunction[interaction_,pack_,compiledFunction_] := Block[{res,inp},
	inp = (Input/.pack)/.(Args/.interaction);
	PrintTemporary["compiled function input = "<>ToString[inp//FullForm] ];
	res = compiledFunction@@inp;
	PrintTemporary[Output/.pack];
    PrintTemporary[res];
];


testCaseCompare[interaction_,pack_]:=Module[{cf},
	cf = testCaseCreateCompiledFunction[interaction,pack];
	PrintTemporary["compiled function="<>ToString[cf]];
	testCaseCompareCompiledFunction[interaction,pack,cf]
];


ToExpression[testCaseExpressionsClear];


(* ::Section:: *)
(*Compare the results generated by plugging the trial values into the Mathematica generated expressions to the results from the list of terms*)


matchExpEval[exp_,eval_] := Module[{i,expNorm,evalNorm,dev},
	expNorm = Table[Symbol[ToString[exp[[i]][[2]]]]->exp[[i]][[1]],{i,1,Length[exp]}];
	evalNorm=Table[Symbol[ToString[eval[[i]][[2]]]]->eval[[i]][[1]],{i,1,Length[eval]}];
	If[(expNorm===evalNorm),
		(* then *) PrintTemporary["They match!"];,
		(* else *)PrintTemporary["The match was not perfect"];
		For[i=1,i<=Length[expNorm],i++,
			dev = expNorm[[i]][[2]]-evalNorm[[i]][[2]];
			If[dev!=0,
				(*then*) PrintTemporary[ToString[expNorm[[i]][[1]]]<>" deviation = "<>ToString[dev]];
			];
		];
	];
];


(* ::Text:: *)
(*Generate a tree representation of the expression list*)


initializeTree[] := Module[{},
Clear[tree,treeSpan];
tree[i_] ={Height->0,TreeSpan->0,Branch->666,Description->0,SubBranches->{},Color->Black};
treeSpan = Table[1,{200}];
treeValues = {};
treeHeight = 0;
];


locateBranches[term_] := Module[
	{poss,tl,i},
	poss={};
	tl = term;
	If[Length[tl]==0,tl={tl}];
	For[i=1,i<=Length[tl],i++,
		AppendTo[poss,tree[tl[[i]]]];
	];
	Return[poss];
];


createBranch[height_,branch_,descr_,subBranches_,color_] := Module[{},
	(*PrintTemporary["createBranch::"];*)
	tree[branch] = {
		Height->height, 
		TreeSpan->treeSpan[[height]], 
		Branch->branch,
		Description->descr, 
		SubBranches->subBranches,
		Color->color
	};
	treeSpan[[height]] = treeSpan[[height]]+1;
	(*PrintTemporary["createBranch:: treeHeight=",treeHeight];*)
	treeHeight = Max[treeHeight,height];
	(*PrintTemporary["treeHeight=",treeHeight];*)
	AppendTo[treeValues,branch];
	Return[tree[branch]];
];


locateCreateBranches[term_] := Module[{i,mx,loc,branch,branches},
loc = locateBranches[term];
(*PrintTemporary["loc=",loc];*)
mx = Max[Table[Height/.loc[[i]],{i,1,Length[loc]}],2];
branches={};
For[i=1,i<=Length[loc],i++,
(*PrintTemporary["locateCreateBranches:: loc[[i]] =", loc[[i]]];*)
tbranch = (Branch/.loc[[i]]);
(*PrintTemporary["locateCreateBranches:: testing if tbranch==666, tbranch=",tbranch];*)
branch = loc[[i]];
If[tbranch==666,
(*PrintTemporary["locateCreateBranches:: tbranch is 666"];*)
branch=createBranch[mx-1,term[[i]],ToString[term[[i]]],{},Black];
];
(*PrintTemporary["locateCreateBranches:: branch=",branch];*)
AppendTo[branches,branch];
];
Return[branches];
];


createRoot[rule_] := Module[
	{lhs,rhs,tl,list,syms,branches,mx,branchNames,i,col,head},
	lhs = rule[[1]];
	rhs = rule[[2]];
	If[Length[lhs]==0,tl={lhs},tl=Table[lhs[[i]],{i,Length[lhs]}];];
	syms=Cases[tl,_Symbol];
	nums=Complement[tl,syms];
	(*PrintTemporary["term=", term[[1]]," syms=",syms," nums=",nums];*)
	branches = locateCreateBranches[syms];
	(*PrintTemporary["branches=",branches];*)
	mx = Max[Table[Height/.branches[[i]],{i,1,Length[branches]}]];
	branchNames = Table[Branch/.branches[[i]],{i,1,Length[branches]}];
	head = Head[lhs];
	desc = ToString[rhs]<>":"<>ToString[head];
	col = Red;
	If[head==Times,col=Green];
	If[head==Plus,col=Blue];
	If[head==Symbol,col=Black];
	If[head==power2,col=Green];
	If[Length[nums]>0,desc=desc<>"("<>ToString[nums[[1]]//InputForm]<>")";];
	createBranch[mx+1,rhs,desc,branchNames,col];
];


createEntireTree[allTerms_]:=Module[{i,terms},
initializeTree[];
terms = selectRules[allTerms];
Map[createRoot,terms];
];


heightOfTree[] := treeHeight


treePos[side_,el_,asx_,asy_] := Module[{max,dev,x,y},
max = Max[treeSpan];
dev = max/treeSpan[[Height/.el]];
x = ((Height/.el)-1);
y = dev*(((TreeSpan/.el)-1)+0.5);
Return[{(side*0.3+x)*asx,y*asy}];
];


branchesInto[x_] := Module[{refs,i,br},
refs={};
For[i=1,i<=Length[treeValues],i++,
br = tree[treeValues[[i]]];
If[MemberQ[SubBranches/.br,x],AppendTo[refs,Branch/.br]];
];
Return[refs];
];


branchesFrom[x_] := SubBranches/.tree[x]


buildTreeGraphics[asx_,asy_]:=Module[
	{ha,subBranches,tv,tvv,br,i,j,lines,graphics,color,line,txt},
	graphics = {};
	lines={};
	ha=0.1;
	For[i=1,i<=Length[treeValues],i++,
		tv=treeValues[[i]];
		br = tree[tv];
		color = Color/.br;
	    If[color==Color,color=Purple];
		subBranches = SubBranches/.br;
		For[j=1,j<=Length[subBranches],j++,
			tvv = tree[subBranches[[j]]];
			line = {color,Line[{treePos[-1,br,asx,asy],treePos[1,tvv,asx,asy]}]};
			AppendTo[lines,line];
		];
		txt = {color,Text[Description/.br,treePos[0,br,asx,asy],Background->Automatic]};
		AppendTo[graphics,txt];
	];
	Return[Flatten[{lines,graphics}]];
];


graphRules[rules_] := Module[{gr},
	createEntireTree[rules];
	gr = buildTreeGraphics[0.1,10];
	Show[Graphics[gr],AspectRatio->0.41]
];


packGraph[pack_] := graphRules[packSelectRules[pack]]


(* ::Section:: *)
(*Simplification Rules*)


simplifyRules = { 
Power[uu_,2]->power2[uu], 
Power[uu_,2.]->power2[uu], 
Power[uu_,3]->power2[uu]*uu,
Power[uu_,4]->power2[uu]*power2[uu],
Power[uu_,Rational[1,2]]->mysqrt[uu],
(*Power[uu_,Rational[-1,2]]->mysqrt[reciprocal[uu]], (*Slower, results in more square roots*) *)
Power[uu_,Rational[-1,2]]->reciprocal[mysqrt[uu]], (*Faster, results in fewer square roots *)
Power[uu_,-1.5]->1.0/(uu*Sqrt[uu]),
Power[uu_,-2.5`]->1.0/(uu*uu*Sqrt[uu]),
Power[uu_,-3.5`]->1.0/(uu*uu*uu*Sqrt[uu]),
Power[uu_,Rational[-3,2]]->1.0/(power1[uu]*mysqrt[uu]),
Power[uu_,Rational[-5,2]]->1.0/(power2[uu]*mysqrt[uu]),
Power[uu_,Rational[-7,2]]->1.0/(power2[uu]*power1[uu]*mysqrt[uu]),
Power[uu_, 1.]->uu,
power1[uu_]->uu,
Power[uu_,-1]->reciprocal[uu],
Power[uu_,-6]->(powern3[uu]*powern3[uu]),
Power[uu_,-2]->(powern2[uu]),
Power[uu_,-3]->(powern2[uu]*reciprocal[uu]),
Power[uu_,-7]->(powern4[uu]*powern3[uu]),
Power[uu_,-8]->(powern4[uu]*powern4[uu]),
Power[uu_,-4]->(powern2[uu]*powern2[uu]),
Power[uu_,-5]->(powern3[uu]*powern2[uu]),
powern2[uu_] -> power2[reciprocal[uu]],
powern4[uu_] -> power2[powern2[uu]],
powern3[uu_] -> powern2[uu]*reciprocal[uu],
Times[1.`,u1_,u2_] -> Times[u1,u2],
Times[1.`,u1_] -> u1,
Times[2,u1_]->2.0*u1,
Times[-2,u1_]->-2.0*u1,
Times[2.,u1_] ->2.0*u1,
Times[-2.,u1_]->-2.0*u1,
Times[6,u1_]->6.0*u1,
Times[-6,u1_]->-6.0*u1,
Times[6.,u1_] ->6.0*u1,
Times[-6.,u1_]->-6.0*u1,
Times[12,u1_]->12.0*u1,
Times[-12,u1_]->-12.0*u1,
Times[12.,u1_] ->12.0*u1,
Times[-12.,u1_]->-12.0*u1,
Times[Rational[-1,2],uu_]->-0.5*uu,
Times[Rational[1,2],uu_]->0.5*uu
};


findRuleForSymbolInRuleList[sym_, ruleList_] := Select[ruleList, Head[#] == Rule && #[[2]] == sym &];
packFindRuleForSymbol[rulePack_, sym_] := Select[packSelectRules[rulePack], #[[2]] == sym &];


(* ::Subsection:: *)
(*Now combine the optimizations into a single function*)


TimesOptimize = True;
PlusOptimize = True;


packOptimize[pack_] := Module[
	{np,outs,rules,code,result,mathExp,fileNamePrefix, fileName},
	uniqueVariableIndex = 1;
	outs = Output/.pack;
	np = pack;
	Print[Style["Set TimesSimplify and PlusSimplify to turn these simplifications off and on",Blue,Italic,24]];
	Print[Style["PlusOptimize = "<>ToString[PlusOptimize],Blue,Italic,20]];
	Print[Style["TimesOptimize = "<>ToString[TimesOptimize],Blue,Italic,20]];
	Print["Collecting terms"]
	np = collectTerms[np,outs]/.simplifyRules;
	np = collectTerms[np,outs]/.simplifyRules;
	np = collectTerms[np,outs]/.simplifyRules;
	np = collectTerms[np,outs]/.simplifyRules;
	If[TimesOptimize==True,
		np = timesSimplify[np,outs]/.simplifyRules;
		,
		PrintTemporary[Style["Skipping timesSimplify",Blue,Italic,20]];
	];
	np = collectTerms[np,outs]/.simplifyRules;
	If[PlusOptimize==True,
		np = plusSimplify[np,outs]/.simplifyRules;
		,
		PrintTemporary[Style["Skipping plusSimplify",Blue,Italic,20]];
	];	
	rules = Rules/.np;
	cFormWriteToFiles[rules,pack];
(*	code = cFormListOfTerms[rules,outs,(AdditionalCDeclares/.pack)];

	fileName = fileNamePrefix<>".cc";
	PrintTemporary["Writing C-code to file("<>fileName<>")"];
	st = OpenWrite[fileName];
	WriteString[st,code];
	Close[st];
'*)
	fileNamePrefix = "_"<>(Name/.pack)<>"_math";
	mathExp = packAsMathematicaBlock[rules,outs];
	result = {
				Name->(Name/.pack),
				AdditionalCDeclares->(AdditionalCDeclares/.pack),
				Input->(Input/.pack),
				Output->outs,
				EnergyFunction->(EnergyFunction/.pack),
				DerivativeVariables->(DerivativeVariables/.pack),
				Rules->(Rules/.np),
				HessianStructure->(HessianStructure/.pack),
				CCode->"NOT USED",
				MathCode->mathExp
			};
	PrintTemporary["Mathematica expression"];
	PrintTemporary[result];
	PrintTemporary["Writing mathematica expression to file:"];
	PrintTemporary["File: "<>fileNamePrefix];
	PrintTemporary["Directory:"<>ToString[Directory[]]];
	Put[result,fileNamePrefix<>".m"];
	Return[result];
];


toArgs[str_] := StringReplace[str,{"{"->"(","}"->")"}];


energyFunctionName[pack_]:=Block[{},
	Return["_evaluateEnergyOnly_"<>(Name/.pack)];
];


writeFiniteDifferenceForce[sout_,pack_] := Block[
	{inputs,derivs,i,vars1,vars2,svars1,svars2,d,str},
	inputs = Input/.pack;
	derivs = DerivativeVariables/.pack;
	For[i=1,i<=Length[derivs],i++,
		d = derivs[[i]];
		vars1 = inputs/.{d->d-delta2};
		vars2 = inputs/.{d->d+delta2};
		svars1 = toArgs[ToString[vars1]];
		svars2 = toArgs[ToString[vars2]];
		str="TEST_FORCE( ";
		str = str<>energyFunctionName[pack];
		str = str<>", ";
		str = str<>" delta2*2.0,\n";
		str = str<>svars1;
		str = str<>",\n";
		str = str<>svars2;
		str = str<>",\n";
		str = str<>oneForceName[d];
		str = str<>",\n";
		str = str<>"index";
		str = str<>");\n";
		WriteString[sout,str];
	];
];


assembleFiniteDifferenceDiagonalHessian[pack_,d_] := Block[
	{res,inputs,derivs,vars1,vars2,vars3,svars1,svars2,svars3},
	inputs = Input/.pack;
	derivs = DerivativeVariables/.pack;
	vars1 = inputs/.{d->d-delta2};
	vars2 = inputs;
	vars3 = inputs/.{d->d+delta2};
	svars1 = toArgs[ToString[vars1]];
	svars2 = toArgs[ToString[vars2]];
	svars3 = toArgs[ToString[vars3]];
	res="TEST_DIAGONAL_HESSIAN( ";
	res=res<>energyFunctionName[pack]<>", ";
	res=res		<>"delta2*2.0,\n";
	res=res<>svars1<>",\n";
	res=res<>svars2<>",\n";
	res=res<>svars3<>",\n";
	res=res<>oneDiagonalHessianName[d]<>",\n";
	res=res<>"index";
	res=res<>");\n";
	Return[res];
];


assembleFiniteDifferenceOffDiagonalHessian[pack_,d1_,d2_] := Block[
	{res,inputs,derivs,vmm,vmp,vpm,vpp,smm,smp,spm,spp},
	inputs = Input/.pack;
	derivs = DerivativeVariables/.pack;
	vmm = inputs/.{d1->d1-delta2,d2->d2-delta2};
	vmp = inputs/.{d1->d1-delta2,d2->d2+delta2};
	vpm = inputs/.{d1->d1+delta2,d2->d2-delta2};
	vpp = inputs/.{d1->d1+delta2,d2->d2+delta2};
	smm = toArgs[ToString[vmm]];
	smp = toArgs[ToString[vmp]];
	spm = toArgs[ToString[vpm]];
	spp = toArgs[ToString[vpp]];
	res="TEST_OFF_DIAGONAL_HESSIAN( "
			<>energyFunctionName[pack]<>",\n"
			<>"delta2*2.0,\n"
			<>smm<>",\n"
			<>spm<>",\n"
			<>smp<>",\n"
			<>spp<>",\n"
			<>oneOffDiagonalHessianName[ToString[d1]<>ToString[d2]]<>",\n"
			<>"index"
			<>");\n";
	Return[res];
];


writeFiniteDifferenceHessian[sout_,pack_] := Block[
	{inputs,derivs,i1,i2,d1,d2,tst},
	inputs = Input/.pack;
	derivs = DerivativeVariables/.pack;
	For[i1=1,i1<=Length[derivs],i1++,
		For[i2=i1,i2<=Length[derivs],i2++,
			d1 = derivs[[i1]];
			d2 = derivs[[i2]];
			If[i1==i2,
				tst=assembleFiniteDifferenceDiagonalHessian[pack,d1];
				,
				tst=assembleFiniteDifferenceOffDiagonalHessian[pack,d1,d2];
				,
				tst = "UNKNOWN";
			];
			WriteString[sout,tst];
		];
	];
];


writeFiniteDifferenceDebug[sout_,pack_]:= Module[
	{},
	WriteString[sout,"// FINITE_DIFFERENCE_TEST \n"];
	writeFiniteDifferenceForce[sout,pack];
	writeFiniteDifferenceHessian[sout,pack];
];
SetAttributes[writeFiniteDifferenceDebug,HoldFirst];


Clear[writeSerializeCode]


writeSerializeCode[pack_, vars_] := Module[{fileName,st,var,i},
	fileName = "_"<>(Symbol["Name"]/.pack)<>"_debugEvalSerialize.cc";
	st = OpenWrite[fileName];
	Print["Writing xml output debug code to: ", fileName ];
	WriteString[st,"// XML_EVAL\n"];
	For [ i=1, i<=Length[vars], i++,
		var = vars[[i]];
		WriteString[st,"\tnode->attributeIfNotDefault(\""<>
						"eval"<>ToString[var]<>"\", this->eval."<>ToString[var]<>
						", 0.0 );\n "];
	];
	Close[st];
];


writeOutputVariablesForDebugging[pack_]:=Block[
	{var,vars,inps,st,i,outs,fileName},

	inps = Input/.pack;
	outs = Output/.pack;
    vars = Flatten[{inps,outs}];
	fileName = "_"<>(Name/.pack)<>"_debugFiniteDifference.cc";
	st = OpenWrite[fileName];
	Print["Writing finite difference debug code to: ", fileName ];
	writeFiniteDifferenceDebug[st,pack];
    Close[st];
	fileName = "_"<>(Name/.pack)<>"_debugEvalDeclares.cc";
	st = OpenWrite[fileName];
	Print["Writing debug variable declares to: ", fileName ];
	WriteString[st,"// DEFINE_EVAL\n"];
	WriteString[st,"\tstruct {\n"];
	For [ i=1, i<=Length[vars], i++,
		var = vars[[i]];
		WriteString[st,"\tdouble\t"<>ToString[var]<>";\n"];
	];
	WriteString[st,"\t} eval;\n"];
	Close[st];

	writeSerializeCode[pack,vars];

	fileName = "_"<>(Name/.pack)<>"_debugEvalSet.cc";
	st = OpenWrite[fileName];
	Print["Writing set variables debug code to: ", fileName ];
	WriteString[st,"// SET_EVAL \n "];
	For [ i=1, i<=Length[vars], i++,
		var = vars[[i]];
		WriteString[st,"\tEVAL_SET("<>
						ToString[var]<>", "<>ToString[var]<>
						");\n "];
	];
	Close[st];
];
	


(* ::Subtitle:: *)
(*Helper functions*)


Assign[var_,exp_]:=Rule[exp,var];


(* Extract the variable name information into a list of rules *)
ExtractNames[names_] := {
	Name->names[[1]],
	XYZ->names[[2]],
	Index->names[[3]],
	Offset->names[[4]],
	Func->If[Length[names]>=5,names[[5]],0]
};


EnergyAccumulate[macroPrefix_,energyName_]:=CCode[macroPrefix<>"_ENERGY_ACCUMULATE("<>ToString[energyName]<>");"];


ForceAccumulate[macroPrefix_,nameRules_,forceName_]:= CCode[
	macroPrefix<>"_FORCE_ACCUMULATE("<>ToString[Index/.nameRules]
	<>", "<>ToString[Offset/.nameRules]<>", "<>ToString[forceName]<>" );"];


AppendGradientAndForce[macroPrefix_,be_,outputs_,rawEnergyFn_,rawNames_] := Module[
	{i,gradName,zforceName,energyFn,nameInfo,varName,names,deriv},
	PrintTemporary["About to append Gradient and Force terms"];
	names = Evaluate[rawNames];
	energyFn = Evaluate[rawEnergyFn];
	For[i=1,i<=Length[names],i++,
		nameInfo = ExtractNames[names[[i]]];
		varName = Name/.nameInfo;
		gradName = Symbol["g"<>ToString[varName]];
		zforceName = Symbol["f"<>ToString[varName]];
		deriv = D[energyFn,varName];
		PrintTemporary["Appending force term: "<>ToString[zforceName]];
		AppendTo[be,Assign[gradName,deriv]];
		AppendTo[be,Assign[zforceName,-gradName]];
		AppendTo[outputs,zforceName];
		AppendTo[be,ForceAccumulate[macroPrefix,nameInfo,zforceName]];
	];
];
SetAttributes[AppendGradientAndForce,HoldAll];


HessianAccumulate[macroPrefix_,name1Info_,name2Info_,hessianName_]:= CCode[
	macroPrefix<>"_HESSIAN_ACCUMULATE("
	<>ToString[Index/.name1Info]<>", "
	<>ToString[Offset/.name1Info]<>", "
	<>ToString[Index/.name2Info]<>", "
	<>ToString[Offset/.name2Info]<>", "
	<>ToString[hessianName]<>");"];


AppendOneHessianElement[macroPrefix_,be_,outputs_,hessianStructure_,energyFn_,names_,i_,j_]:=Module[
	{name1Info,name2Info,var1Name,var2Name,hessianName,hess,prefix,zfirstDeriv},
	name1Info = ExtractNames[names[[i]]];
	name2Info = ExtractNames[names[[j]]];
	var1Name = Name/.name1Info;
	var2Name = Name/.name2Info;
	If[i==j,prefix="d"(*Diagonal*),prefix="o"];
	hessianName = Symbol[prefix<>"h"<>ToString[var1Name]<>ToString[var2Name]];
	zfirstDeriv = D[energyFn,var1Name];
	hess = D[zfirstDeriv,var2Name];
	AppendTo[be,Assign[hessianName,hess]];
	AppendTo[outputs,hessianName];
	AppendTo[be,HessianAccumulate[macroPrefix,name1Info,name2Info,hessianName]];
	PrintTemporary["Appended hessian element ("<>ToString[j]
				<>","<>ToString[i]<>") for variable: "<>ToString[hessianName]];
	hessianStructure[[i,j]]=Length[outputs];
	If[i!=j,hessianStructure[[j,i]]=Length[outputs]];
];
SetAttributes[AppendOneHessianElement,HoldAll];


AppendHessianDiagonalElements[macroPrefix_,be_,outputs_,hessianStructure_,rawEnergyFn_,rawNames_] := Module[
	{names,energyFn,i},
	names = Evaluate[rawNames];
	energyFn = Evaluate[rawEnergyFn];
	For[i=1,i<=Length[names],i++,
		AppendOneHessianElement[macroPrefix<>"_DIAGONAL",be,outputs,hessianStructure,energyFn,names,i,i];
	];
];
SetAttributes[AppendHessianDiagonalElements,HoldAll];


AppendHessianOffDiagonalElements[macroPrefix_,be_,outputs_,hessianStructure_,rawEnergyFn_,rawNames_] := Module[
	{energyFn,name1Info,name2Info,var1Name,var2Name,hessianName,col,row,names},
	names = Evaluate[rawNames];
	energyFn = Evaluate[rawEnergyFn];
	For[col=1,col<Length[names],col++,
		For[row=col+1,row<=Length[names],row++,
			AppendOneHessianElement[macroPrefix<>"_OFF_DIAGONAL",be,outputs,hessianStructure,energyFn,names,col,row];
		];
	];
];
SetAttributes[AppendHessianOffDiagonalElements,HoldAll];


(* Now do everything *)
AppendGradientForceAndHessian[macroPrefix_,be_,outputs_,hessianStructure_,rawEnergyFn_,rawVarNames_]:= Module[
	{energyFn,varNames},
	energyFn = Evaluate[rawEnergyFn];
	varNames = Evaluate[rawVarNames];
	AppendTo[be,CCode["#ifdef "<>macroPrefix<>"_CALC_FORCE //["]];
	AppendTo[be,CCode["if ( calcForce ) {"] ];
	AppendGradientAndForce[macroPrefix,be,outputs,energyFn,varNames];
	AppendTo[be,CCode["#ifdef "<>macroPrefix<>"_CALC_DIAGONAL_HESSIAN //["]];
	AppendTo[be,CCode["if ( calcDiagonalHessian ) {"]];
	AppendHessianDiagonalElements[macroPrefix,be,outputs,hessianStructure,energyFn,varNames];
	AppendTo[be,CCode["#ifdef "<>macroPrefix<>"_CALC_OFF_DIAGONAL_HESSIAN //["]];
	AppendTo[be, CCode["if ( calcOffDiagonalHessian ) {" ]];
	AppendHessianOffDiagonalElements[macroPrefix, be, outputs,hessianStructure, energyFn, varNames];
	AppendTo[be, CCode["} /*if calcOffDiagonalHessian */ "]];
	AppendTo[be, CCode["#endif /* "<>macroPrefix<>"_CALC_OFF_DIAGONAL_HESSIAN ]*/"]];
	AppendTo[be,CCode["} /*calcDiagonalHessian */"]];
	AppendTo[be, CCode["#endif /* "<>macroPrefix<>"_CALC_DIAGONAL_HESSIAN ]*/"]];
	AppendTo[be,CCode["} /*calcForce */"]];
	AppendTo[be, CCode["#endif /* "<>macroPrefix<>"_CALC_FORCE ]*/"]];
];
SetAttributes[AppendGradientForceAndHessian,HoldAll];


(* For rigid body calculations *)
AppendGradientForceIfCalcForce[macroPrefix_,be_,outputs_,rawEnergyFn_,rawVarNames_]:= Module[
	{energyFn,varNames},
	energyFn = Evaluate[rawEnergyFn];
	varNames = Evaluate[rawVarNames];
	AppendTo[be,CCode["#ifdef "<>macroPrefix<>"_CALC_FORCE //["]];
	AppendTo[be,CCode["if ( calcForce ) {"] ];
	AppendGradientAndForce[macroPrefix,be,outputs,energyFn,varNames];
	AppendTo[be,CCode["} /*calcForce */"]];
	AppendTo[be, CCode["#endif /* "<>macroPrefix<>"_CALC_FORCE ]*/"]];
];
SetAttributes[AppendGradientForceIfCalcForce,HoldAll];


EndPackage[]



