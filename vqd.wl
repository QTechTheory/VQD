(* ::Package:: *)

BeginPackage["VQD`"];


Needs["QuEST`"];


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
TrappedIonHub::usage="Returns device specification of a multi-nodes Trapped ions based on the device built by the QCSHub.";
TrappedIonInnsbruck::usage="Returns device specification of a string of Trapped ions base on the device built by the University of Innsbruck.";

(*rydberg quantum devices/neutral atoms.*)
RydbergHub::usage="Returns device specification of a Rydberg/Neutral Atom device based on the device built by the QCSHub.";
RydbergWisconsin::usage="Returns device specification of a Rydberg/Neutral Atom device based on the device built by the University of Wisconsin.";

(*nuclear-vacancy center devices.*)
NVCenterDelft::usage="Returns device specification of a Nitrogen-Vacancy diamond center device based on the device built by the University of Delft.";
NVCenterHub::usage="Returns device specification of a Nitrogen-Vacancy diamond center device based on the device built by the QCSHub.";


ParameterDevices::usage="Show all parameters used to specify all devices. To see parameters related to a device, e.g., NVCenterHub, use Options[NVCenterHub].";


BeginPackage["`ParameterDevices`"];


(* Parameters *)
BField::usage="The electromagnetic field strength in the z-direction from the lab reference with unit Tesla.";
DurTwoGate::usage="Duration of two qubit gates with rotation of \[Pi]";
DurMeas::usage="Duration of measurement";
DurInit::usage="Duration of initialisation";
DurShuffle::usage="Duration to shuffle location of ions";
entangling::usage="Crosstalk error model pre equation (4) on applying the XX M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
ErrCT::usage="Error coefficient of the crosstalk with entanglement model";
ErrSS::usage="Error coefficient of the crosstalk with stark shift model";
EFSingleXY::usage="Error fraction/ratio, {depolarising, dephasing} of the single qubit X and Y rotations. Sum of the ratio must be 1 or 0 (off).";
EFCZ::usage="Error fraction/ratio, {depolarising, dephasing} of the controlled-Z gates. Sum of the ratio must be 1 or 0 (off).";
EFCRotXY::usage="Error fraction/ratio, {depolarising, dephasing} of the controlled-Rx and -Ry gates. Sum of the ratio must be 1 or 0 (off).";
ExchangeRotOn::usage="Maximum interaction j on the passive qubit crosstalk when applying CZ gates; The noise form is C[Rz[j.\[Theta]]] It must be a square matrix with size (nqubit-2)x(nqubit-2).";
ExchangeRotOff::usage="Crosstalks error C-Rz[ex] on the passive qubits when not applying two-qubit gates.";
FidSingleXY::usage="Fidelity(ies) of single Rx[\[Theta]] and Ry[\[Theta]] rotations obtained by random benchmarking.";
FidCRotXY::usage="The fidelity of controlled-X and controlled-Y gates.";
FidMeas::usage="Fidelity of measurement";
FidInit::usage="Fidelity of qubit initialisation";
FidCZ::usage="Fidelity(ies) of the CZ gates.";
FreqSingleXY::usage="Rabi frequency(ies) for the single X- and Y- rotations with unit MHz";
FreqCRotXY::usage="Rabi frequency for the controlled-X and controlled-Y rotations with unit MHZ.";
FreqCZ::usage="Rabi frequency(ies) for the CZ gate with unit MHz.";
MSCrossTalk::usage="(entangling OR starkshift) The crosstalk model in applying M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
NIons::usage="The total number of ions in a trapped ion device.";
qubitsNum::usage="The number of physical active qubits for computations.";
RabiFreq::usage="The Rabi frequency frequency in average or on each qubit with unit MHz.";
starkshift::usage="Crosstalk error model using stark shift on  applying the Exp[-i\[Theta]XX], namely M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
T1::usage="T1 duration(s) in \[Mu]s. Exponential decay time for the state to be complete mixed.";
T2::usage="T2 duration(s) in \[Mu]s. Exponential decay time for the state to be classical with echo applied.";
T2s::usage="T2* duration(s) in \[Mu]s. Exponential decay time for the state to be classical.";


BField::error="`1`";
ExchangeRotOff::error="`1`";
ExchangeRotOn::error="`1`";
FidCRotXY::error="`1`";
FidCZ::error="`1`";
FidSingleXY::error="`1`";
FidCRotXY::warning="`1`";
FidCZ::warning="`1`";
FidSingleXY::warning="`1`";

FreqCRotXY::error="`1`";
FreqCZ::error="`1`";
Larmor::error="`1`";
qubitsNum::error="`1`";
T1::error="`1`";
T2::error="`1`";
T2s::error="`1`";


(* Custom gates *)
SWAPLoc::usage="Swap the spatial locations of two qubits";
Wait::usage="Wait gate, doing nothing";
Init::usage="Initialise qubit to state |0>";


(*Visualisations*)
DrawIons::usage="Draw the current string of ions";


EndPackage[];


(*All definitions of modules*)


Begin["`Private`"];


validate::usage="Validate expression, throw error if false.";
validate[value_,expr_,err_,msg_,format_:nothing]:=(If[expr[format[value]],format[value],Throw[Message[err::error,msg]]])

nothing[arg___]:=arg
checkAss::usage="check if it's an association with length len.";
checkAss[ass_,len_]:=AssociationQ[ass]&&Length[ass]===len&&And@@NumberQ/@Values@ass
checkAss[ass_,len_,f_]:=AssociationQ[ass]&&Length[ass]===len&&And@@f/@Values@ass
num2Ass[arg_Real,len_Integer]:=<|Table[i->arg,{i,0,-1+len}]|>
num2Ass[arg_Integer,len_Integer]:=<|Table[i->arg,{i,0,-1+len}]|>
num2Ass[arg_Association,len_Integer]:=arg


numass[len_]:="not a number or association of numbers with length "<>ToString[len]
fidass[len_]:="not a fidelity number or association of fidelities with length "<>ToString[len]


(*noise forms on the SiliconDelft*)
exCZON::usage="exCZON[target_qubit,qubitNum,exchangeMatrix]. Exchange C-Rz[j] interaction when CZ gate on";
exCZON[targ_,nq_,ex_]:=Subscript[C, #-1][Subscript[Rz, #][ex[[targ,#]]]]&/@Delete[Range[nq-1],targ]


SiliconDelft[OptionsPattern[]]:=Module[
{
(*validate and format parameter specification*)
(*Numbers*)
qubitsnum=Catch@validate[OptionValue@qubitsNum,IntegerQ,qubitsNum,"not an integer"],
bfield=Catch@validate[OptionValue@BField,NumberQ,BField,"not a number."],
(*Fractions*)
efsinglexy=Catch@validate[OptionValue@EFSingleXY,(Total[#]==1||Total[#]==0)&,EFSingleXY,"not a fraction with total 1 "],
efcz=Catch@validate[OptionValue@EFCZ,(Total[#]==1||Total[#]==0)&,EFCZ,"not a fraction with total 1 "],
(*Number as average or association to specify each*)
t1=Catch@validate[OptionValue@T1,checkAss[#,OptionValue@qubitsNum]&,T1,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
t2=Catch@validate[OptionValue@T2,checkAss[#,OptionValue@qubitsNum]&,T2,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
t2s=Catch@validate[OptionValue@T2s,checkAss[#,OptionValue@qubitsNum]&,T2s,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
rabifreq=Catch@validate[OptionValue@RabiFreq,checkAss[#,OptionValue@qubitsNum]&,Larmor,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
exchangerotoff=Catch@validate[OptionValue@ExchangeRotOff,checkAss[#,-1+OptionValue@qubitsNum]&,ExchangeRotOff,numass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
freqcz=Catch@validate[OptionValue@FreqCZ,checkAss[#,-1+OptionValue@qubitsNum]&,FreqCZ,numass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
freqcrotxy=Catch@validate[OptionValue@FreqCRotXY,checkAss[#,-1+OptionValue@qubitsNum]&,FreqCRotXY,numass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
freqsinglexy=Catch@validate[OptionValue@FreqSingleXY,checkAss[#,OptionValue@qubitsNum]&,FreqSingleXY,numass[OptionValue@qubitsNum],num2Ass[#,OptionValue@qubitsNum]&],
(*Number as average or association to specify each fidelity*)
fidcrotxy=Catch@validate[OptionValue@FidCRotXY,checkAss[#,-1+OptionValue@qubitsNum,0<=#<=1&]&,FidCRotXY,fidass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
fidcz=Catch@validate[OptionValue@FidCZ,checkAss[#,-1+OptionValue@qubitsNum,0<=#<=1&]&,FidCZ,fidass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
fidsinglexy=Catch@validate[OptionValue@FidSingleXY,checkAss[#,OptionValue@qubitsNum,0<=#<=1&]&,FidSingleXY,fidass[OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
(*a matrix*)
exchangeroton=Catch@validate[OptionValue@ExchangeRotOn,SquareMatrixQ,ExchangeRotOn,StringForm["not a square matrix with dim ``^2",-1+OptionValue@qubitsNum]],
(*probability error of depolarising and dephasing noise *)
er1xy=fid2DepolDeph[#,OptionValue@EFSingleXY,1,FidSingleXY,True]&/@OptionValue[FidSingleXY],
ercz=fid2DepolDeph[#,OptionValue@EFCZ,2,FidCZ,True]&/@OptionValue[FidCZ],
er2xy=fid2DepolDeph[#,OptionValue@EFCRotXY,2,FidCRotXY,True]&/@OptionValue[FidCRotXY],

(* frequently used stuff *)
qubits=Range[0,-1+OptionValue@qubitsNum]
},

Module[
{\[CapitalDelta]t},
<|
(*no hidden qubits/ancilla here *)
DeviceDescription -> "Delft Silicon device with "<>ToString[qubitsnum]<>"-qubits arranged as a linear array with nearest-neighbor connectivity and control qubits are the lower ones.",
NumAccessibleQubits -> qubitsnum,
NumTotalQubits -> qubitsnum,
Aliases -> {
	Subscript[Wait, q__][t_] :> {}
	},	
Gates ->{
	Subscript[Wait, -1][t_]:><|
		NoisyForm->{},
		GateDuration->t
	|>,
(* Singles *)
	Subscript[Rz, q_][\[Theta]_]:> <|
		NoisyForm->{Subscript[Rz, q][\[Theta]]},
		GateDuration->0(*virtual gate, instant, noiseless*)
		|>,
	Subscript[Rx,q_][\[Theta]_]:><|
		NoisyForm->{Subscript[Rx, q][\[Theta]],Subscript[Depol, q][er1xy[q][[1]]],Subscript[Deph, q][er1xy[q][[2]]],Sequence@@Table[Subscript[U, j][OffResRabiOsc[rabifreq[q],rabifreq[j],Abs[\[Theta]]]],{j,Delete[qubits,q+1]}]},
		GateDuration->Abs[\[Theta]]/(2\[Pi] freqsinglexy[q]) 
	|>,
	Subscript[Ry,q_][\[Theta]_]:><|
		NoisyForm->{Subscript[Ry, q][\[Theta]],Subscript[Depol, q][er1xy[q][[1]]],Subscript[Deph, q][er1xy[q][[2]]],Sequence@@Table[Subscript[U, j][OffResRabiOsc[rabifreq[q],rabifreq[j],Abs[\[Theta]]]],{j,Delete[qubits,q+1]}]},
		GateDuration->Abs[\[Theta]]/(2\[Pi] freqsinglexy[q])
	|>,
(* Twos *)
	Subscript[C, p_][Subscript[Z, q_]]/; (q-p)===1  :><|
		(*The last bit undo the exchange in the passive noise *)
		NoisyForm->{Subscript[C, p][Subscript[Z, q]],Subscript[Depol, p,q][ercz[p][[1]]],Subscript[Deph, p,q][ercz[p][[2]]],Sequence@@exCZON[q,qubitsnum,exchangeroton],Sequence@@Table[Subscript[C, i][Subscript[Rz, i+1][(-\[CapitalDelta]t/\[Pi])*exchangerotoff[i]]],{i,Complement[qubits,{p,q,qubitsnum-1}]}]}, 
		GateDuration->0.5/freqcz[p] 
	|>,
		Subscript[C, p_][Subscript[Rx, q_][\[Theta]_]]/; (q-p)===1  :><|
		NoisyForm->{Subscript[C, p][Subscript[Rx, q][\[Theta]]],Subscript[Depol, p,q][ercz[p][[1]]],Subscript[Deph, p,q][ercz[p][[2]]]}, 
		GateDuration->0.5/freqcz[p] 
	|>
	
	
},
(* Declare that \[CapitalDelta]t will refer to the duration of the current gate/channel. *)
DurationSymbol -> \[CapitalDelta]t, 
(* Passive noise *)
Qubits :> {
		q_/;(q==qubitsnum-1) :> <|
		PassiveNoise ->{Subscript[Deph, q][.5(1-E^(-\[CapitalDelta]t/t2s[q]))],Subscript[Depol, q][.75(1-E^(-\[CapitalDelta]t/t1[q]))]}
		|>
		,
		q_/;(q<qubitsnum-1) :> <|
		PassiveNoise ->{Subscript[C, q][Subscript[Rz, q+1][(\[CapitalDelta]t/\[Pi])*exchangerotoff[q]]],Subscript[Deph, q][.5(1-E^(-\[CapitalDelta]t/t2s[q]))],Subscript[Depol, q][.75(1-E^(-\[CapitalDelta]t/t1[q]))]}
		|>
				
		}
|>
]
]


OffResRabiOsc::usage="OffResRabiOsc[rabi_freq_active, rabi_freq_passive, duration] Off resonant Rabi Oscillation.";
OffResRabiOsc[\[CapitalOmega]_,\[CapitalOmega]0_,t_]:=Module[{\[CapitalDelta],\[CapitalOmega]R},
\[CapitalDelta]=\[CapitalOmega]-\[CapitalOmega]0;
\[CapitalOmega]R=Sqrt[\[CapitalOmega]^2+\[CapitalDelta]^2];
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


End[];


EndPackage[];


Needs["VQD`ParameterDevices`"]
