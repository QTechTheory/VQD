(* ::Package:: *)

BeginPackage["VQD`"];

(* some constants *)
(*Planck constant, Js*)
hbar=1.054571817*10^34 ;
(*Bohr magneton J/T*)
\[Mu]B=9.274009994*10^-24;

(* Devices *)
(*silicon devices*)
SiliconDelft::usage="Returns device specification of a Silicon device based on the device built by the University of Delft.";
SiliconHub::usage="Returns devices specification of a Silicon device based on the device built by the QCSHub.";

(*superconducting qubit devices*)
SuperconductingFZJ::usage="Returns device specification of a Superconducting qubit device based on the device built by Forschungzentrum Juelich.";
SuperconductingHub::usage="Returns device specification of a Superconducting qubit device based on the device built by the QCSHub.";

(*trapped ion devices*)
TrappedIonOxford::usage="Returns device specification of a multi-nodes Trapped ions based on the device built by the Oxford/Hub.";
TrappedIonInnsbruck::usage="Returns device specification of a string of Trapped ions base on the device built by the University of Innsbruck.";

(*rydberg quantum devices/neutral atoms.*)
RydbergHub::usage="Returns device specification of a Rydberg/Neutral Atom device based on the device built by the QCSHub.";
RydbergWisconsin::usage="Returns device specification of a Rydberg/Neutral Atom device based on the device built by the University of Wisconsin.";

(*nuclear-vacancy center devices.*)
NVCenterDelft::usage="Returns device specification of a Nitrogen-Vacancy diamond center device based on the device built by the University of Delft.";
NVCenterHub::usage="Returns device specification of a Nitrogen-Vacancy diamond center device based on the device built by the QCSHub.";

ParameterDevices::usage="Show all parameters used to specify all devices. To see parameters related to a device, e.g., NVCenterHub, use Options[NVCenterHub].";

(* Other functions *)
PartialTrace::usage="PartialTrace[qureg, qubits_to_be_traced_out]. Return the partial trace as a matrix.";
RandomMixState::usage="RandomMixState[nqubits, nsamples:None]. Return a random mixed quantum density state matrix. It is obtained by the sum random unitary matrices.
By default, the random unitaries sample is \!\(\*SuperscriptBox[\(2\), \(\(nqubit\)\(\\\ \)\)]\)";
SerializeCircuit::usage="SerializeCircuit[circuit]. Execute the gates without parallelism with exception. For instance, in the case of neutral atoms, gates parallelism can be done when the spatial distance is outside blockade radius. This function should be executed before InsertCiruitNoise.";
(* Custom gates *)
SWAPLoc::usage="Swap the spatial locations of two qubits";
Wait::usage="Wait gate, doing nothing";
Init::usage="Initialise qubit to state |0>";
CZ::usage="Controlled-Z operation";
Ent::usage="Remote entanglement operation";
SplitZ::usage="SplitZ[node, zone_destination]. Split a string of ions in a zone of a trapped-ion Oxford device";
CombS::usage="CombS[node, zone_destination]. Combine a string of ions to a zone of a trapped-ion Oxford device";

(*Visualisations*)
DrawIons::usage="Draw the current string of ions";
BeginPackage["`ParameterDevices`"];
(* Parameters *)
BField::usage="The electromagnetic field strength in the z-direction from the lab reference with unit Tesla.";
DurTwoGate::usage="Duration of two qubit gates with rotation of \[Pi]";
DurMeas::usage="Duration of measurement";
DurInit::usage="Duration of initialisation";
DurRead::usage="Readout duration in \[Mu]s";
DurShuffle::usage="Duration to shuffle location of ions";
DurMove::usage="Duration for physically moving operation in Trapped Ions such as SplitZ and CombS.";
DDActive::usage="Apply dynamical decoupling: use T2 in the model if set True, otherwise use T2* if set False.";
entangling::usage="Crosstalk error model pre equation (4) on applying the XX M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
ErrCT::usage="Error coefficient of the crosstalk with entanglement model";
ErrSS::usage="Error coefficient of the crosstalk with stark shift model";
EFSingleXY::usage="Error fraction/ratio, {depolarising, dephasing} of the single qubit X and Y rotations. Sum of the ratio must be 1 or 0 (off).";
EFCZ::usage="Error fraction/ratio, {depolarising, dephasing} of the controlled-Z gates. Sum of the ratio must be 1 or 0 (off).";
EFInit::usage="Error fraction/ratio {depolarising, dephasing} of the initialisation gate.  Sum of the ratio must be 1 or 0 (off).";
EFCRot::usage="Error fraction/ratio, {depolarising, dephasing} of the controlled-Rx and -Ry gates. Sum of the ratio must be 1 or 0 (off).";
EFEnt::usage="Error fraction/ratio {depolarising, dephasing} of remote entanglement.  Sum of the ratio must be 1 or 0 (off).";
EFRead::usage="Error fraction/ratio of {depolarising,dephasing} of the readout. Sum of the ratio must be 1 or 0 (off)." ;
ExchangeRotOn::usage="Maximum interaction j on the passive qubit crosstalk when applying CZ gates; The noise form is C[Rz[j.\[Theta]]] It must be a square matrix with size (nqubit-2)x(nqubit-2).";
ExchangeRotOff::usage="Crosstalks error C-Rz[ex] on the passive qubits when not applying two-qubit gates.";
FidSingleXY::usage="Fidelity(ies) of single Rx[\[Theta]] and Ry[\[Theta]] rotations obtained by random benchmarking.";
FidEnt::usage="Fidelity of remote entanglement operation.";
FidMeas::usage="Fidelity of measurement";
FidInit::usage="Fidelity of qubit initialisation";
FidCZ::usage="Fidelity(ies) of the CZ gates.";
FreqSingleXY::usage="Rabi frequency(ies) for the single X- and Y- rotations with unit MHz";
FreqCZ::usage="Rabi frequency(ies) for the CZ gate with unit MHz.";
FreqEnt::usage="Frequency of remote entanglement.";
FidRead::usage="Readout fidelity";
MSCrossTalk::usage="(entangling OR starkshift) The crosstalk model in applying M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
Nodes::usage="Entire nodes of a trapped ions system <|node1 -> number_of_qubits_1, ... |>";
NIons::usage="The total number of ions in a trapped ion device.";
OffResonantRabi::usage="Put the noise due to off-resonant Rabi oscillation when applying single qubit rotations.";
ParallelGates::usage="Gates that are executed in parallel even with the call SerializeCircuit[]";
qubitsNum::usage="The number of physical active qubits for computations.";
QubitFreq::usage="The Qubit frequency for each qubit with unit MHz.";
RabiFreq::usage="The Rabi frequency frequency in average or on each qubit with unit MHz.";
RepeatRead::usage="The number of repeated readout (n) performed. The final fidelity of readout is \!\(\*SuperscriptBox[\(FidRead\), \(n\)]\).";
Meas::usage="Perform measurement on the qubits";
starkshift::usage="Crosstalk error model using stark shift on  applying the Exp[-i\[Theta]XX], namely M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
ShowNodes::usage="Draw all Ions on every nodes within the zones";
T1::usage="T1 duration(s) in \[Mu]s. Exponential decay time for the state to be complete mixed.";
T2::usage="T2 duration(s) in \[Mu]s. Exponential decay time for the state to be classical with echo applied.";
T2s::usage="T2* duration(s) in \[Mu]s. Exponential decay time for the state to be classical.";
StdPassiveNoise::usage="Set to True/False. Use the standard passive noise that involves T1, T2 or T2s inputs.";
BField::error="`1`";
DurInit::error="`1`";
DDActive::error="`1`";
ExchangeRotOff::error="`1`";
ExchangeRotOn::error="`1`";
EFRead::error="`1`";
EFInit::error="`1`";
FidCRotXY::error="`1`";
FidCZ::error="`1`";
FidSingleXY::error="`1`";
FidMeas::error="`1`";
FidInit::error="`1`";
FreqCZ::error="`1`";
QubitFreq::error="`1`";
Nodes::error="`1`";
qubitsNum::error="`1`";
OffResonantRabi::error="`1`";
RabiFreq::error="`1`";
StdPassiveNoise::error="`1`";
T1::error="`1`";
T2::error="`1`";
T2s::error="`1`";
FidCRotXY::warning="`1`";
FidCZ::warning="`1`";
FidSingleXY::warning="`1`";
FidMeas::warning="`1`";
FidInit::warning="`1`";
EndPackage[];

(*******All definitions of modules****)
Begin["`Private`"];

Needs["QuEST`"];
Needs["QuEST`Option`"];
Needs["QuEST`Gate`"];
Needs["QuEST`DeviceSpec`"];

numass[len_]:="not a number or association of numbers with length "<>ToString[len]
fidass[len_]:="not a fidelity number or association of fidelities with length "<>ToString[len]
nothing[arg___]:=arg
validate::usage="Validate expression, throw error if false.";
validate[value_,expr_,err_,msg_,format_:nothing]:=(If[expr[format[value]],format[value],Throw[Message[err::error,msg]]])
checkAss::usage="check if it's an association with length len.";
checkAss[ass_,len_]:=AssociationQ[ass]&&Length[ass]===len&&And@@NumberQ/@Values@ass
checkAss[ass_,len_,f_]:=AssociationQ[ass]&&Length[ass]===len&&And@@f/@Values@ass
num2Ass[arg_Real,len_Integer]:=<|Table[i->arg,{i,0,-1+len}]|>
num2Ass[arg_Integer,len_Integer]:=<|Table[i->arg,{i,0,-1+len}]|>
num2Ass[arg_Association,len_Integer]:=arg

OffResRabiOsc::usage="OffResRabiOsc[rabi_freq, qubit_freq_diff, duration] Off resonant Rabi Oscillation.";
OffResRabiOsc[\[CapitalOmega]_,\[CapitalDelta]_,t_]:=With[{\[CapitalOmega]R=Sqrt[\[CapitalOmega]^2+\[CapitalDelta]^2]},
E^(I \[CapitalDelta] t/2) {{Cos[\[CapitalOmega]R*t/2]-I*Sin[\[CapitalOmega]R*t/2]*\[CapitalDelta]/\[CapitalOmega]R,  -I*Sin[\[CapitalOmega]R*t/2]*\[CapitalOmega]/\[CapitalOmega]R            },
			  {-I*Sin[\[CapitalOmega]R*t/2]*\[CapitalOmega]/\[CapitalOmega]R,             Cos[\[CapitalOmega]R*t/2]+I Sin[\[CapitalOmega]R*t/2]*\[CapitalDelta]/\[CapitalOmega]R}}				 					  			  
]

fid2DepolDeph::usage = "fid2DepolDeph[totFid, {ratio.Depol,ratio.Deph}, nqerr, err, avgfid:True]. 
Return the parameter for depolarizing and dephasing noise that gives the total fidelity totFid [0,1], where totFid is the average fidelity obtained from random benchmarking (set avgfid:False if it's the worst fidelity)
correction is a constant adjusted to params of dephasing to get average fidelity from the worst fidelity.
";
fid2DepolDeph[totfid_, errratio_, nqerr_, errval_, avgfid_:True] := Module[
  {sol, rate, pdepol, pdeph, entfid},
  If[
  (*set to perfect fidelity of zero error *)
  (totfid==1||0==Total@errratio),
  {pdepol,pdeph}={0,0}
  ,
    
  If[
  (** 1-qubit error **)
  nqerr === 1, 
  sol=If[avgfid,
		(*estimate parameters from entanglement fidelity*)
		entfid=(3totfid-1)/2;
		Solve[1-rate*errratio[[1]]-rate*errratio[[2]]+4/3 errratio[[1]]*errratio[[2]]*rate^2==entfid,{rate}]
		,
		(*estimate parameters from the worst fidelity: from |+\[RightAngleBracket] state*)
		Solve[1-2errratio[[1]]*rate/3-errratio[[2]]*rate+4*errratio[[1]]*errratio[[2]]*rate^2/3==totfid,{rate}]
	];
	(* check validity of numbers, set it to the worst parameter if exceeds *)
	{pdepol, pdeph}= errratio*Min@Abs[rate/.sol];
	If[pdepol>0.75,
		Message[errval::warning,StringForm["(warning) fidelity is too low; 1-qubit depolarization parameter is ``. Set it to 3/4.",pdepol]];
		pdepol=0.75;
		];
	If[pdeph>0.5,
		Message[errval::warning,StringForm["(warning) fidelity is too low; 1-qubit dephasing parameter is ``. Set it to 1/2.",pdeph]];
		pdeph=0.5;
	];
	,
	(** 2-qubits error **)
	sol=If[avgfid,
	(*estimate parameters from entanglement fidelity*)
	entfid=(5*totfid-1)/4;
	Solve[1-rate*errratio[[1]]-rate*errratio[[2]]+16/15*rate^2*errratio[[1]]*errratio[[2]]==entfid,{rate}]
	,
	(*estimate parameters from the worst fidelity: from |++\[RightAngleBracket] state*)
	Solve[1-errratio[[2]]*rate-errratio[[1]]*rate*4/5+errratio[[1]]*errratio[[2]]*rate^2*16/15==totfid,{rate}]
	];
	(* check validity of numbers, set it to the worst fidelity if exceeds *)
	{pdepol, pdeph}= errratio*Min@Abs[rate/.sol];
	If[pdepol>15/16,
		Message[errval::warning,StringForm["(warning) Fidelity might be too low; 2-qubit depolarization parameter is ``. Set it to 15/16.",pdepol]];
		pdepol=15/16;
	];
		
	If[pdeph>3/4,
		Message[errval::warning,StringForm["(warning) Fidelity might be too low; 2-qubit dephasing parameter is ``. Set it to 3/4.",pdeph]];
		pdeph=3/4;
	];
];
];	
  {pdepol, pdeph}
]
grouptwo::usage="grouptwo[list], group a list into two elements";
grouptwo[list_]:=ReplaceList[Sort@list,{p___,a_,b_,q___}:>{a,b}]
(*additional error models *)
bitFlip::usage="Return kraus operator with bitflip error";
bitFlip[fid_,q_]:=With[{e=1-fid},
	Subscript[Kraus, q][{Sqrt[1-e]*{{1,0},{0,1}},Sqrt[e]*{{1,0},{0,1}}}]
	]

bitFlip[fid_,p_,q_]:=With[{e=1-fid},
	Subscript[Kraus, p,q][{
		Sqrt[1-e]*IdentityMatrix[4],
		Sqrt[e/3]*{{0,1,0,0},{1,0,0,0},{0,0,0,1},{0,0,1,0}},
		Sqrt[e/3]*{{0,0,1,0},{0,0,0,1},{1,0,0,0},{0,1,0,0}},
		Sqrt[e/3]*{{0,0,0,1},{0,0,1,0},{0,1,0,0},{1,0,0,0}}
}]]

(***** SILICON_DELFT *****)
SiliconDelft[OptionsPattern[]]:=With[
{

(*validate and format parameter specification*)
(*Numbers*)
qubitsnum=Catch@validate[OptionValue@qubitsNum,IntegerQ,qubitsNum,"not an integer"],
(*Fractions*)
efsinglexy=Catch@validate[OptionValue@EFSingleXY,(Total[#]==1||Total[#]==0)&,EFSingleXY,"not a fraction with total 1 "],
efcz=Catch@validate[OptionValue@EFCZ,(Total[#]==1||Total[#]==0)&,EFCZ,"not a fraction with total 1 "],

(*Number as average or association to specify each*)
t1=Catch@validate[OptionValue@T1,checkAss[#,OptionValue@qubitsNum]&,T1,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
t2=Catch@validate[OptionValue@T2,checkAss[#,OptionValue@qubitsNum]&,T2,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
t2s=Catch@validate[OptionValue@T2s,checkAss[#,OptionValue@qubitsNum]&,T2s,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
rabifreq=Catch@validate[OptionValue@RabiFreq,checkAss[#,OptionValue@qubitsNum]&,RabiFreq,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
qubitfreq=Catch@validate[OptionValue@QubitFreq,checkAss[#,OptionValue@qubitsNum]&,QubitFreq,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
freqcz=Catch@validate[OptionValue@FreqCZ,checkAss[#,-1+OptionValue@qubitsNum]&,FreqCZ,numass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],

(*Number as average or association to specify each fidelity*)
fidcz=Catch@validate[OptionValue@FidCZ,checkAss[#,-1+OptionValue@qubitsNum,0<=#<=1&]&,FidCZ,fidass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
fidsinglexy=Catch@validate[OptionValue@FidSingleXY,checkAss[#,OptionValue@qubitsNum,0<=#<=1&]&,FidSingleXY,fidass[OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],

(*single things*)
fidread=Catch@validate[OptionValue@FidRead,0<=#<=1&,FidRead,"invalid fidelity"],
durread=Catch@validate[OptionValue@DurRead,NumberQ,DurMeas,"invalid duration"],
repeatread=Catch@validate[OptionValue@RepeatRead,IntegerQ,RepeatRead,"not an integer"],

(* assoc or boolean *)
exchangerotoff=Catch@validate[OptionValue@ExchangeRotOff,Or[AssociationQ@#,#===False]&,ExchangeRotOff,"Set to association or False"],

(* True/False*)
ddactive=Catch@validate[OptionValue@DDActive,BooleanQ,DDActive,"Set to true or false"],
offresonantrabi=Catch@validate[OptionValue@OffResonantRabi,BooleanQ,OffResonantRabi,"Set to true or false"],
stdpassivenoise=Catch@validate[OptionValue@StdPassiveNoise,BooleanQ,StdPassiveNoise,"Set to true or false"],


(*a matrix or a boolean*)
exchangeroton=Catch@validate[OptionValue@ExchangeRotOn,Or[SquareMatrixQ[#],BooleanQ[#]]&,ExchangeRotOn,StringForm["set to false or specify with a square matrix with dim ``^2",-1+OptionValue@qubitsNum]],

(*probability error of depolarising and dephasing noise *)
er1xy=fid2DepolDeph[#,OptionValue@EFSingleXY,1,FidSingleXY,True]&/@OptionValue[FidSingleXY],
ercz=fid2DepolDeph[#,OptionValue@EFCZ,2,FidCZ,True]&/@OptionValue[FidCZ],

(* frequently used stuff *)
qubits=Range[0,-1+OptionValue@qubitsNum]
},

Module[
{\[CapitalDelta]t, miseq,initf,measf, passivenoisecirc, offresrabi, stdpn, exczon,durinit,sroterr},

stdpn[q_,dd:True]:=If[stdpassivenoise,{Subscript[Deph, q][.5(1-E^(-\[CapitalDelta]t/t2[q]))],Subscript[Depol, q][.75(1-E^(-\[CapitalDelta]t/t1[q]))]},{}];
stdpn[q_,dd:False]:=If[stdpassivenoise,{Subscript[Deph, q][.5(1-E^(-\[CapitalDelta]t/t2s[q]))],Subscript[Depol, q][.75(1-E^(-\[CapitalDelta]t/t1[q]))]},{}];
passivenoisecirc[q_Integer]:=Flatten@{If[q<qubitsnum-1 && AssociationQ@exchangerotoff,Subscript[C, q][Subscript[Rz, q+1][(\[CapitalDelta]t/\[Pi])*exchangerotoff[q]]],{}],stdpn[q,stdpassivenoise]};

(*single rotation noise *)
offresrabi[q_,\[Theta]_]:=If[offresonantrabi,Table[Subscript[U, j][OffResRabiOsc[rabifreq[q],qubitfreq[j]-qubitfreq[q],Abs[\[Theta]]]],{j,Delete[qubits,q+1]}],{}];

(*Exchange rotation C-Rz[j] interaction when CZ gate on*)
exczon[targ_]:=If[ListQ@exchangeroton,Subscript[C, #-1][Subscript[Rz, #][exchangeroton[[targ,#]]]]&/@Delete[Range[qubitsnum-1],targ],{}];

(*Errors on single rotations*)
sroterr[q_,\[Theta]_]:=Flatten@{Subscript[Depol, q][er1xy[q][[1]]],Subscript[Deph, q][er1xy[q][[2]]],offresrabi[q,\[Theta]]};
(*measurement and initialisation sequence*)
miseq[q__]:=(Length@{q}>1)&&((Sort[{q}]===Range[0,Max@q])||(Sort[{q}]===Range[Min@q,-1+qubitsnum]));
(* final init state is 1000...0001: 2 reads+1 cond-X *)
		
initf[q__]:=Which[MemberQ[{q},0],(*start*)
		Flatten@{
		(*mimics 2 readout *)
		{Subscript[Damp, 0][1-(1-fidread)^2],Subscript[X, 0],Subscript[Damp, 1][1-(1-fidread)^2]},
		(* perfect partial swaps + 1 readout *)
		Table[{Subscript[C, i][Subscript[X, i-1]],sroterr[i-1,\[Pi]],Subscript[C, i-1][Subscript[X, i]],sroterr[i,\[Pi]]},{i,Complement[{q},{0,1}]}]
		,
		{Subscript[Damp, 1][1-(1-fidread)]}
		}
		,
		MemberQ[{q},qubitsnum-1],
		(*end*)
		Flatten@{
		{Subscript[Damp, qubitsnum-1][1-(1-fidread)^2],Subscript[X, qubitsnum-1],Subscript[Damp, qubitsnum-2][1-(1-fidread)^2]},
		(* perfect partial swaps + 1 readout *)
		Table[{Subscript[C, i][Subscript[X, i+1]],sroterr[i+1,\[Pi]],Subscript[C, i+1][Subscript[X, i]],sroterr[i,\[Pi]]},{i,Complement[{q},{qubitsnum-1,qubitsnum-2}]}]
		,
		{Subscript[Damp, qubitsnum-2][1-(1-fidread)]}
		}
		,
		True,
		Throw[Message[Init::error,"Error in Init operation"]]
];

		
(*duration of a single initialisation: 3 readout + 2X + 1 CX *)	
durinit[q__]:=durread*3+1.5*Total@@Table[rabifreq[i],{i, Complement[{q},{0,qubitsnum-1}]}]; 
						
(* measurment has bitflip errors at the edge and single qubit errors in the middle *)		
measf[q__]:=Which[
	MemberQ[{q},0],(*start*)
	Flatten@{
	(*mimics 2 readout  *)
	{bitFlip[1-(1-fidread)^2,0,1],Subscript[M, 0],Subscript[M, 1],Subscript[Damp, 0][1-(1-fidread)^2],Subscript[X, 0],Subscript[Damp, 1][1-(1-fidread)^2]},
	(* 2 rotation errors + 1 readout *)
	Table[{Subscript[Depol, i][er1xy[i][[1]]],Subscript[Deph, i][er1xy[i][[2]]],Subscript[M, i],Subscript[Damp, i][1-(1-fidread)]},{i,Complement[{q},{0,1}]}]
	}
	,
	MemberQ[{q},qubitsnum-1],(*end*)
	Flatten@{
	(*mimics 2 readout  *)
	{bitFlip[1-(1-fidread)^2,qubitsnum-1,qubitsnum-2],Subscript[M, qubitsnum-1],Subscript[M, qubitsnum-2],Subscript[Damp, qubitsnum-1][1-(1-fidread)^2],Subscript[X, qubitsnum-1],Subscript[Damp, qubitsnum-2][1-(1-fidread)^2]},
	(* 2 rotation errors + 1 readout *)
	Table[{Subscript[Depol, i][er1xy[i][[1]]],Subscript[Deph, i][er1xy[i][[2]]],Subscript[M, i],Subscript[Damp, i][1-(1-fidread)]},{i,Complement[{q},{qubitsnum-1,qubitsnum-2}]}]
	}
	,
	True,
	Throw[Message[M::error,"Error in Measurement operation"]]
];

<|
(*no hidden qubits/ancilla here *)
DeviceDescription -> "Delft Silicon device with "<>ToString[qubitsnum]<>"-qubits arranged as a linear array with nearest-neighbor connectivity and control qubits are the lower ones.",
NumAccessibleQubits -> qubitsnum,
NumTotalQubits -> qubitsnum,

Aliases -> {
	Subscript[Wait, q__][t_] :> {},
	Subscript[Init, q__]:> {},
	Subscript[Meas, q__]:> {}
	}
	,	
Gates ->{
	Subscript[Wait, q__][t_]:><|
		NoisyForm->Flatten@Table[passivenoisecirc[i],{i,Flatten@{q}}],
		GateDuration->t
	|>,
(* Measurements and initialisation *)
	Subscript[Meas, q__]/; miseq[q] :><|
		NoisyForm-> measf[q],
		GateDuration->durread
	|>,
	Subscript[Init, q__]/; miseq[q] :><|
		NoisyForm-> initf[q],
		GateDuration-> durinit[q]
	|>,		
(* Singles *)
	Subscript[Rx,q_][\[Theta]_]:><|
		NoisyForm->Flatten@{Subscript[Rx, q][\[Theta]],sroterr[q,\[Theta]]},
		GateDuration->Abs[\[Theta]]/(2\[Pi] rabifreq[q]) 
	|>,
	Subscript[Ry,q_][\[Theta]_]:><|
		NoisyForm->Flatten@{Subscript[Ry, q][\[Theta]],sroterr[q,\[Theta]]},
		GateDuration->Abs[\[Theta]]/(2\[Pi] rabifreq[q])
	|>,
(* Twos *)
	Subscript[C, p_][Subscript[Z, q_]]/; (q-p)===1  :><|
		(*The last bit undo the exchange in the passive noise *)
		NoisyForm->{Subscript[C, p][Subscript[Z, q]],Subscript[Depol, p,q][ercz[p][[1]]],Subscript[Deph, p,q][ercz[p][[2]]],Sequence@@exczon[q]}, 
		GateDuration->0.5/freqcz[p] 
	|>		
},
(* Declare that \[CapitalDelta]t will refer to the duration of the current gate/channel. *)
DurationSymbol -> \[CapitalDelta]t, 
(* Passive noise *)
Qubits :> {
		q_/;(0<=q<qubitsnum) :> <|
		PassiveNoise ->passivenoisecirc[q]
		|>		 
		}	
	|>
]
]
(***** ENDOF SILICON_DELFT *****)

(***** TRAPPED_IONS_OXFORD *****)
createNodes::usage="createNodes[ <| zone1 -> nq1, zone2 -> nq2,...|>]. Return nodes with 4 zones, its map to qubits register, and the total number of qubits. Initially, all qubits are in zone 1";
createNodes[args__Association]:=Module[{q,qglob=0,nodes,qmap},
nodes=Map[<|1->Range[#],2->{},3->{},4->{}|>&,args];
qmap=Map[<|#[1]/.{q_Integer:>(q->qglob++)}|>&,nodes];
{nodes, qmap,qglob}
]

(** drawing-related trapped ions **)
showIons[nodes_]:=Module[
{qulocs,i,colors},
colors=rainbowcol[Length@nodes];
qulocs=<|Table[node-><|Table[i=0;zone-><|Table[qubit->{i++,0},{qubit,nodes[node][zone]}]|>,{zone,Keys@nodes[node]}]|>,{node,Keys@nodes}]|>;
i=0;
Grid[Table[
i++;{node,Sequence@@Table[drawZone[qulocs[node][zone],"zone "<>ToString@zone,colors[[i]]],{zone,Keys@nodes[node]}]},{node,Keys@nodes}],Spacings->{0,1}]
]
rainbowcol[n_]:=Reverse@Table[ColorData["Rainbow",(i-1)/n],{i,n}]
drawZone[locszone_,label_,color_]:=Show[Graphics[Table[{Thick,color,Sphere[loc,.25]},{loc,Values@locszone}]],Graphics[Table[Text[k,locszone[k]],{k,Keys@locszone}]],
ImageSize->{50*Max[1,Length@locszone],50}, Frame->True,FrameTicks->False,FrameStyle->Black,FrameLabel->label,FrameMargins->Tiny]

(* zone-related functions *)
getZone[q_,node_]:=(Association@Flatten@Table[Table[v->k,{v,node[k]}],{k,Keys@node}])[q]

(* Legitimate split move *)
legSplit[nodes_,nodename_,q1_,q2_,zone_]:=With[{z1=getZone[q1,nodes[nodename]],z2=getZone[q2,nodes[nodename]],node=nodes[nodename]},
If[z1!=z2,Return@False];
If[Abs[Position[node[z1],q1][[1,1]]- Position[node[z1],q2]][[1,1]]!=1,Return@False];
If[Length@Flatten[Table[node[z],{z,Range[Min[z1,zone]+1,Max[z1,zone]-1]}]]>0,Return@False];
True
]

(*Legitimate combine moves *)
legComb[nodes_,nodename_,q1_,q2_,zone_]:=Module[{node=nodes[nodename],z1,z2,cond1,cond2,zstart,ps,pz,qs,qz,cond3},
	(* combine to the one of the zone *)
	z1=getZone[q1,node];z2=getZone[q2,node];
	cond1=Or[z1===zone,z2===zone];
	If[\[Not]cond1,Return[False]];
	If[z1===zone,
		zstart=z2;qs=q2;qz=q1,
		zstart=z1;qs=q1;qz=q2];
	(* merge from start zone to zone destination *)
	ps=Position[node[zstart],qs][[1,1]];
	pz=Position[node[zone],qz][[1,1]];
	cond2=Which[
	(* move down to a higher zone *)
	zstart<zone,
		Length@node[zstart][[ps+1;;]]+Length@node[zone][[;;pz-1]],
	zstart>zone,
	(*move up to a lower zone *)
		Length@node[zone][[pz+1;;]]+Length@node[zstart][[;;ps-1]],
	zstart===zone,
	(* doesn't move *)
	Length@node[zone][[Min[ps,pz]+1;;Max[ps,pz]-1]]
	];
	If[cond2!=0,Return[False]];

	(* no qubits in between zones *)
	cond3=If[zone===zstart,
	True,
	Length@Flatten[Table[node[z],{z,Range[Min[zstart,zone]+1,Max[zstart,zone]-1]}]]===0];
	If[\[Not]cond3,Return[False]];
	True
]

(* Legitimate split move *)
Subscript[splitZ, i_,j_][inodes_,node_,zone_]:=Module[{zq1,zq2,pos,sq,szone,nodes},
	nodes=inodes;
	zq1=getZone[i,nodes@node];
	zq2=getZone[j,nodes@node];
	Assert[zq1===zq2,"both qubits must have the same zone"];
	pos=<|i->Position[nodes[node][zq1],i][[1,1]],j->Position[nodes[node][zq1],j][[1,1]]|>;
	sq=Keys@pos; (*sorted qubits by its post: {low,high}*)

	(*split the zone*)
	szone=TakeDrop[nodes[node][zq1],Values[pos][[1]]];
	If[zone>zq1,
		(* move to a higher zone, lower qubit remains *)
		nodes[node][zq1]=szone[[1]];
		nodes[node][zone]=Join[szone[[2]],nodes[node][zone]],
		(* move to a lower zone, high qubit remains *)
		nodes[node][zone]=Join[nodes[node][zone],szone[[1]]];
		nodes[node][zq1]=szone[[2]];
	];
	nodes
]

Subscript[CombS, i_,j_][inodes_,nodename_,zone_]:=Module[{nodes,node,z1,z2,ps,pz,sq,zstart,qs,qz},
	nodes=inodes;
	node=nodes[nodename];
	z1=getZone[i,node];
	z2=getZone[j,node];
	If[z1===zone,
		zstart=z2;qs=j;qz=i,
		 zstart=z1;qs=i;qz=j
	];
	(* merge from start zone to zone destination *)
	ps=Position[node[zstart],qs][[1,1]];
	pz=Position[node[zone],qz][[1,1]];
	Which[
		(* move down to a higher zone *)
	zstart<zone,
		nodes[nodename][zone]=Join[{qs},node[zone]];
		nodes[nodename][zstart]=node[zstart][[;;-2]];
	,
	zstart>zone,
		(*move up to a lower zone *)
		nodes[nodename][zone]=Join[node[zone],{qs}];
		nodes[nodename][zstart]=node[zstart][[2;;]];
	];
	nodes
]

(*
Zone1 :prepare, store, detect
Zone2 :prepare, store, detect, logic
Zone3 :prepare, store, detect, logic
Zone4 :remote entangle
*)
TrappedIonOxford[OptionsPattern[]]:=With[
{
initnodes=OptionValue@Nodes,
durinit=OptionValue@DurInit,
durmove=OptionValue@DurMove,
durread=OptionValue@DurRead,
efent=OptionValue@EFEnt,
fidinit=OptionValue@FidInit,
fidread=OptionValue@FidRead,
freqcz=OptionValue@FreqCZ,
freqent=OptionValue@FreqEnt,
fident=OptionValue@FidEnt,
rabifreq=OptionValue@RabiFreq,
t1=OptionValue@T1,
t2=OptionValue@T2,

(*probability error of depolarising and dephasing noise *)
er1xy=Association[#->fid2DepolDeph[OptionValue[FidSingleXY][#],OptionValue[EFSingleXY][#],1,FidSingleXY,True]&/@Keys[OptionValue@Nodes]], 
ercz=Association[#->fid2DepolDeph[OptionValue[FidCZ][#],OptionValue[EFCZ][#],2,FidCZ,True]&/@Keys[OptionValue@Nodes]],
erent=fid2DepolDeph[OptionValue@FidEnt,OptionValue@EFEnt,2,FidEnt,True]
},
Module[
{\[CapitalDelta]t, nodes, qmap, qnum, passivenoise,checkpsd,checklog,checkrent,gnoise,swaploc},
{nodes,qmap,qnum}=createNodes[initnodes];

(*TODO:get standard passive noise when certain qubits are acted on*)
passivenoise[node_,t_,q__]:=Flatten@Table[{Subscript[Depol, qmap[node][i]][.75(1-E^(-t/t1[node]))],Subscript[Deph, qmap[node][i]][.5(1-E^(-t/t2[node]))]},{i,Complement[Flatten@Values[nodes[node]],{q}]}];
passivenoise[node_,t_]:=Flatten@Table[{Subscript[Depol, qmap[node][i]][.75(1-E^(-t/t1[node]))],Subscript[Deph, qmap[node][i]][.5(1-E^(-t/t2[node]))]},{i,Flatten@Values[nodes[node]]}];

(** check zones **)
(* prepare, store, detect *)
checkpsd[q_,node_]:=If[MemberQ[{1,2,3},getZone[q,nodes[node]]],True,False];
(*logic*)
checklog[q_,node_]:=If[MemberQ[{2,3},getZone[q,nodes[node]]],True,False];
checklog[p_,q_,node_]:=If[MemberQ[{2,3},getZone[q,nodes[node]]]&&(getZone[q,nodes@node]===getZone[p,nodes@node]),True,False];
(*remote entangle*)
checkrent[q1_,node1_,q2_,node2_]:=If[(4===getZone[q1,nodes[node1]])&&(4===getZone[q2,nodes[node2]]),True,False];

(**gate noise **)
gnoise[q_,node_,\[Theta]_]:=Sequence@@{Subscript[Depol, qmap[node][q]][er1xy[node][[1]]],Subscript[Deph, qmap[node][q]][er1xy[node][[2]]]};
gnoise[p_,q_,node_,\[Theta]_]:=Sequence@@{Subscript[Depol, qmap[node][p],qmap[node][q]][ercz[node][[1]]],Subscript[Deph, qmap[node][p],qmap[node][q]][ercz[node][[2]]]};
gnoise[q1_,q2_,node1_,node2_,t_]:=Sequence@@{Subscript[Depol, qmap[node1][q1],qmap[node2][q2]][erent[[1]]],Subscript[Deph, qmap[node1][q1],qmap[node2][q2]][erent[[2]]]};

(** extra operations **)
swaploc[p_,q_,node_]:=Module[{z=getZone[q,nodes[node]],pos,lst},
	lst=nodes[node][z];
	pos=Flatten@{Position[lst,p],Position[lst,q]};
	lst[[pos]]=lst[[Reverse@pos]];
	nodes[node][z]=lst;
	nodes
];

<|
(*no hidden qubits/ancilla here *)
DeviceDescription -> StringForm["Trapped ion device Oxford style with ``. nodes",Length@nodes],
NumAccessibleQubits -> qnum,
NumTotalQubits -> qnum,
Nodes:>nodes,
ShowNodes:>showIons[nodes],
(* Init, Read, Rx, Ry, C[Z], Ent[node1,node2], SWAP, SplitZ, CombS  *)
Aliases -> {
	Subscript[Wait, q__][node_,t_] :> {},
	Subscript[Init, q_][fid_]:> Subscript[Damp, q][fid],
	Subscript[SplitZ, i_,j_]:> Sequence@@{},
	Subscript[CombS, i_,j_]:>Sequence@@{},
	Subscript[Read, q_]:>Subscript[M, q],
	Subscript[Ent, p_,q_]:>Sequence@@{Subscript[H, p],Subscript[C, p][Subscript[X, q]]},
	Subscript[SWAP, i_,j_][t_]:>Sequence@@{}
	}
	,	
Gates ->{
	Subscript[Rx, q_][node_,\[Theta]_]/;checklog[q,node]:><|
		NoisyForm->{Subscript[Rx, qmap[node][q]][\[Theta]],gnoise[q,node,\[Theta]]},
		GateDuration->rabifreq[node]*Abs[\[Theta]]/(2\[Pi])
	|>,
	Subscript[Ry, q_][node_,\[Theta]_]/;checklog[q,node]:><|
		NoisyForm->{Subscript[Ry, qmap[node][q]][\[Theta]],gnoise[q,node,\[Theta]]},
		GateDuration->rabifreq[node]*Abs[\[Theta]]/(2\[Pi])
	|>,
	Subscript[CZ, i_,j_][node_]/;checklog[i,j,node]:><|
		NoisyForm->{Subscript[C, qmap[node][i]][Subscript[Z, qmap[node][j]]],gnoise[i,j,node,\[Pi]]},
		GateDuration->freqcz[node]*0.5
	|>,
	Subscript[Ent, q1_,q2_][node1_,node2_]/;checkrent[q1,node1,q2,node2]:><|
		NoisyForm-> {Subscript[Ent, qmap[node1][q1],qmap[node2][q2]],gnoise[q1,q2,node1,node2,1/freqent]},
		GateDuration->1/freqent
	|>,
	Subscript[Read, q_][node_]/; checkpsd[q,node] :><|
		NoisyForm->{bitFlip[fidread@node,qmap[node][q]],Subscript[Read, qmap[node][q]]},
		GateDuration->durread[node]
	|>,
	Subscript[Init, q_][node_]/; checkpsd[q,node] :><|
		NoisyForm->{Subscript[Init, qmap[node][q]][fidinit[node]]},
		GateDuration->durinit[node]
	|>,
	Subscript[Wait, q__][node_,t_]:><|
		NoisyForm->passivenoise[node,t],
		GateDuration->t
	|>,
	Subscript[SWAP, i_,j_][node_]/;checklog[i,j,node] :><|
		NoisyForm->{Subscript[SWAP, i,j][durmove[node]]},
		GateDuration->durmove[node],
		UpdateVariables->Function[nodes=swaploc[i,j,node]]
	|>,
	Subscript[SplitZ, i_,j_][node_,zone_]/;legSplit[nodes,node,i,j,zone] :><|
		NoisyForm-> Flatten@{Subscript[SplitZ, qmap[node][i],qmap[node][j]]},
		GateDuration->durmove[node],
		UpdateVariables->Function[nodes=Subscript[splitZ, i,j][nodes,node,zone]]
	|>,
	Subscript[CombS, i_,j_][node_,zone_]/;legComb[nodes,node,i,j,zone]:><|
		NoisyForm->Flatten@{Subscript[CombS, qmap[node][i],qmap[node][j]]},
		GateDuration->durmove[node],
		UpdateVariables->Function[nodes=Subscript[CombS, i,j][nodes,node,zone]]
	|>	
}
,
(* Declare that \[CapitalDelta]t will refer to the duration of the current gate/channel. *)
DurationSymbol -> \[CapitalDelta]t, 
(* Passive noise *)
Qubits :> {}
		
	|>
]
]

RandomMixState[nqubits_]:=Module[{size=2^nqubits,gm,um,dm,id},
(* Random states generation: https://iitis.pl/~miszczak/files/papers/miszczak12generating.pdf *)
(* random ginibre matrix *)
gm=RandomReal[NormalDistribution[0,1],{size,size}]+I RandomReal[NormalDistribution[0,1],{size,size}];
(* random unitatry *)
um=RandomVariate[CircularUnitaryMatrixDistribution[size]];
id=IdentityMatrix[size];
dm=(id+um) . gm . ConjugateTranspose[gm] . (id+ConjugateTranspose[um]);
dm/Tr[dm]//Chop
]

GateIndex::usage="GateIndex[gate], returns {the gate, the indices}.";
GateIndex[gate_]:=Module[{p,q,g,base,idx},
{base,idx}=(gate/.{R[t_,Subscript[X, p_] Subscript[X, q_]]-> {XX, {p,q}}, Subscript[g_, q__][_]-> {g, {Sequence@@q}}, Subscript[g_, q__]->{g,{Sequence@@q}}});
{base,Flatten@idx}
];
SetAttributes[SerializeCircuit,HoldAll]
SerializeCircuit[circuit_]:=Module[{circ=circuit,newcircc,circols,idxcol,incol,g1,g2,idx1,idx2},
newcircc={};
While[Length@circ>0 ,
circols=GetCircuitColumns[circ];
(* update equivalent ordering of the circuit *)
circ=Flatten@circols;
(* get indices partitioned wrt circols *)
idxcol=TakeList[Range[Length@circ],Length@#&/@circols][[1]];

(* eliminate non-legitimate gates of the first column *)
AppendTo[newcircc,{}];
incol=<| #->True & /@idxcol |>;

Table[
If[incol[i1],
(* add gate i1 and eliminate the rest *)
AppendTo[newcircc[[-1]],circ[[i1]]];
Table[
If[incol[[i2]],
{g1,idx1}=GateIndex[circ[[i1]]];
{g2,idx2}=GateIndex[circ[[i2]]];
incol[[i2]]=MemberQ[ParallelGates,(g1|g2)]
];
,{i2,Complement[idxcol,{i1}]}];
]
,{i1,idxcol}];
(* update circuit *)
circ=Delete[circ,{#}&/@Keys@Select[incol,#&]];
];
newcircc
]

(*Partial trace on n-qubit
1) reshape to the tensor:ConstantArray[2,2*n]
2) contract
3) reshape to matrix with dim2^mx2^mwhere m=(n-#contract)
*)
PartialTrace[\[Rho]_Integer,qubits__]:=Module[{\[Rho]mat,tmat,nq,pmat,nfin,pairs},
\[Rho]mat=GetQuregMatrix[\[Rho]];
nq=Log2@Length@\[Rho]mat;
(*tensorize*)
tmat=ArrayReshape[\[Rho]mat,ConstantArray[2,2*nq]];

(* contraction pairs, beware the least significant bit convention!! *)
pairs=Reverse@Table[{i,i+nq},{i,Range[nq]}];
pmat=TensorContract[tmat,pairs[[Sequence@@#]]&/@(1+{qubits})];
nfin=nq-Length@{qubits};
ArrayReshape[pmat,{2^nfin,2^nfin}]//Chop
]

PartialTrace[\[Rho]mat_List,qubits__]:=Module[{tmat,nq,pmat,nfin,pairs},
nq=Log2@Length@\[Rho]mat;
(*tensorize*)
tmat=ArrayReshape[\[Rho]mat,ConstantArray[2,2*nq]];
pairs=Table[{i,i+nq},{i,Range[nq]}];
pmat=TensorContract[tmat,pairs[[Sequence@@#]]&/@(1+{qubits})];
nfin=nq-Length@{qubits};
ArrayReshape[pmat,{2^nfin,2^nfin}]//Chop
]

End[];
EndPackage[];
Needs["VQD`ParameterDevices`"]



