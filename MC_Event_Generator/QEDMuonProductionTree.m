(* ::Package:: *)

(* ::Section:: *)
(*Load FeynCalc and FeynArts*)


If[ $FrontEnd === Null,
		$FeynCalcStartupMessages = False;
		Print["Computation of the matrix element squared for muon production in QED at tree level"];
];
If[$Notebooks === False, $FeynCalcStartupMessages = False];
$LoadFeynArts=True;
<<FeynCalc`
$FAVerbose = 0;


(* ::Section:: *)
(*Generate Feynman diagrams*)


topMuonProd = CreateTopologies[0, 2 -> 2];
diagsMuonProd =
		InsertFields[topMuonProd, {F[2, {1}], -F[2, {1}]} -> {F[2,
		{2}], -F[2, {2}]}, InsertionLevel -> {Classes}, Model -> "SM",
		ExcludeParticles -> {S[1], S[2]}];
Paint[diagsMuonProd, ColumnsXRows -> {1, 1}, Numbering -> None,SheetHeader->None,ImageSize->{256,256}];


(* ::Section:: *)
(*Obtain corresponding amplitudes*)


ampMuonProd = FCFAConvert[CreateFeynAmp[diagsMuonProd, Truncated -> False],
IncomingMomenta->{p1,p2},OutgoingMomenta->{k1,k2},UndoChiralSplittings->True,ChangeDimension->4,List->False,SMP->True];


(* ::Section:: *)
(*Unpolarized process  e^- e^+ -> mu^- mu^+ *)


SetMandelstam[s, t, u, p1, p2, -k1, -k2, SMP["m_e"], SMP["m_e"], SMP["m_mu"], SMP["m_mu"]];
sqAmpMuonProd = (ampMuonProd (ComplexConjugate[ampMuonProd]//FCRenameDummyIndices))//
		PropagatorDenominatorExplicit//Contract//FermionSpinSum[#, ExtraFactor -> 1/2^2]&//
		ReplaceAll[#, DiracTrace :> Tr] &//Contract//Simplify


masslessElectronsSqAmpMuonProd = (sqAmpMuonProd /. {SMP["m_e"] -> 0})//Simplify;


masslessElectronsMuonsSqAmpMuonProd = (masslessElectronsSqAmpMuonProd /. {SMP["m_mu"] -> 0})//Simplify


result = masslessElectronsMuonsSqAmpMuonProd/.s->4 ECM^2/.t->-2 ECM^2(1+Cos[theta])/.u->-2 ECM^2(1-Cos[theta])


res = result/.SMP["cos_W"]->Sqrt[cw2]/.SMP["sin_W"]->Sqrt[sw2]/.SMP["m_Z"]->MZ//Simplify


res1 = res/.cw2->1-sw2/.SMP["e"]->Sqrt[4 \[Pi] alpha]//Simplify


res1>>"/home/guojintseng/Desktop/MC_Event_Generator/calc_dsigma.txt"


res2 = res1/.sw2->0.222246/.ECM->45/.alpha->1/137/.MZ->91/.Cos[theta]->x


Plot[res2,{x,-1,1}]
