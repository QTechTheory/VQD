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
TrappedIonHub::usage="Returns device specification of a multi-nodes Trapped ions based on the device built by the QCSHub.";
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

(*Visualisations*)
DrawIons::usage="Draw the current string of ions";
BeginPackage["`ParameterDevices`"];
(* Parameters *)
BField::usage="The electromagnetic field strength in the z-direction from the lab reference with unit Tesla.";
DurTwoGate::usage="Duration of two qubit gates with rotation of \[Pi]";
DurMeas::usage="Duration of measurement";
DurInit::usage="Duration of initialisation";
DurShuffle::usage="Duration to shuffle location of ions";
DDActive::usage="Apply dynamical decoupling: use T2 in the model if set True, otherwise use T2* if set False.";
entangling::usage="Crosstalk error model pre equation (4) on applying the XX M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
ErrCT::usage="Error coefficient of the crosstalk with entanglement model";
ErrSS::usage="Error coefficient of the crosstalk with stark shift model";
EFSingleXY::usage="Error fraction/ratio, {depolarising, dephasing} of the single qubit X and Y rotations. Sum of the ratio must be 1 or 0 (off).";
EFCZ::usage="Error fraction/ratio, {depolarising, dephasing} of the controlled-Z gates. Sum of the ratio must be 1 or 0 (off).";
EFCRot::usage="Error fraction/ratio, {depolarising, dephasing} of the controlled-Rx and -Ry gates. Sum of the ratio must be 1 or 0 (off).";
EFRead::usage="Error fraction/ratio of {depolarising,dephasing} of the readout. Sum of the ratio must be 1 or 0 (off)." ;
ExchangeRotOn::usage="Maximum interaction j on the passive qubit crosstalk when applying CZ gates; The noise form is C[Rz[j.\[Theta]]] It must be a square matrix with size (nqubit-2)x(nqubit-2).";
ExchangeRotOff::usage="Crosstalks error C-Rz[ex] on the passive qubits when not applying two-qubit gates.";
FidSingleXY::usage="Fidelity(ies) of single Rx[\[Theta]] and Ry[\[Theta]] rotations obtained by random benchmarking.";
FidMeas::usage="Fidelity of measurement";
FidInit::usage="Fidelity of qubit initialisation";
FidCZ::usage="Fidelity(ies) of the CZ gates.";
FreqSingleXY::usage="Rabi frequency(ies) for the single X- and Y- rotations with unit MHz";
FreqCZ::usage="Rabi frequency(ies) for the CZ gate with unit MHz.";
MSCrossTalk::usage="(entangling OR starkshift) The crosstalk model in applying M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
NIons::usage="The total number of ions in a trapped ion device.";
OffResonantRabi::usage="Put the noise due to off-resonant Rabi oscillation when applying single qubit rotations.";
ParallelGates::usage="Gates that are executed in parallel even with the call SerializeCircuit[]";
qubitsNum::usage="The number of physical active qubits for computations.";
RabiFreq::usage="The Rabi frequency frequency in average or on each qubit with unit MHz.";
starkshift::usage="Crosstalk error model using stark shift on  applying the Exp[-i\[Theta]XX], namely M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
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
FidCRotXY::error="`1`";
FidCZ::error="`1`";
FidSingleXY::error="`1`";
FidMeas::error="`1`";
FidInit::error="`1`";
FreqCZ::error="`1`";
Larmor::error="`1`";
qubitsNum::error="`1`";
OffResonantRabi::error="`1`";
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
grouptwo::usage="grouptwo[list], group a list into two elements";
grouptwo[list_]:=ReplaceList[Sort@list,{p___,a_,b_,q___}:>{a,b}]

(***** SILICON_DELFT *****)
SiliconDelft[OptionsPattern[]]:=With[
{

(*validate and format parameter specification*)
(*Numbers*)
qubitsnum=Catch@validate[OptionValue@qubitsNum,IntegerQ,qubitsNum,"not an integer"],
bfield=Catch@validate[OptionValue@BField,NumberQ,BField,"not a number."],
(*Fractions*)
efsinglexy=Catch@validate[OptionValue@EFSingleXY,(Total[#]==1||Total[#]==0)&,EFSingleXY,"not a fraction with total 1 "],
efcz=Catch@validate[OptionValue@EFCZ,(Total[#]==1||Total[#]==0)&,EFCZ,"not a fraction with total 1 "],
efread=Catch@validate[OptionValue@EFRead,(Total[#]==1||Total[#]==0)&,EFRead,"not a fraction with total 1 "],

(*Number as average or association to specify each*)
t1=Catch@validate[OptionValue@T1,checkAss[#,OptionValue@qubitsNum]&,T1,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
t2=Catch@validate[OptionValue@T2,checkAss[#,OptionValue@qubitsNum]&,T2,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
t2s=Catch@validate[OptionValue@T2s,checkAss[#,OptionValue@qubitsNum]&,T2s,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
rabifreq=Catch@validate[OptionValue@RabiFreq,checkAss[#,OptionValue@qubitsNum]&,Larmor,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
freqcz=Catch@validate[OptionValue@FreqCZ,checkAss[#,-1+OptionValue@qubitsNum]&,FreqCZ,numass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
freqsinglexy=Catch@validate[OptionValue@FreqSingleXY,checkAss[#,OptionValue@qubitsNum]&,FreqSingleXY,numass[OptionValue@qubitsNum],num2Ass[#,OptionValue@qubitsNum]&],
(*Number as average or association to specify each fidelity*)
fidcz=Catch@validate[OptionValue@FidCZ,checkAss[#,-1+OptionValue@qubitsNum,0<=#<=1&]&,FidCZ,fidass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
fidsinglexy=Catch@validate[OptionValue@FidSingleXY,checkAss[#,OptionValue@qubitsNum,0<=#<=1&]&,FidSingleXY,fidass[OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],

(*just association*)
fidinit=Catch@validate[OptionValue@FidInit,(And@@Table[0<=j<=1,{j,Values@#}])&,FidInit,"invalid fidelities"],
durinit=Catch@validate[OptionValue@DurInit,(And@@Table[j>=0,{j,Values@#}])&,DurInit,"invalid durations"],
fidmeas=Catch@validate[OptionValue@FidMeas,(And@@Table[0<=j<=1,{j,Values@#}])&,FidMeas,"invalid fidelities"],
durmeas=Catch@validate[OptionValue@DurMeas,(And@@Table[j>=0,{j,Values@#}])&,DurMeas,"invalid durations"],

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
erread=fid2DepolDeph[#,OptionValue@EFRead,2,FidMeas,True]&/@OptionValue[FidMeas],

(* frequently used stuff *)
qubits=Range[0,-1+OptionValue@qubitsNum]
},

Module[
{\[CapitalDelta]t, miseq,initf,measf, passivenoisecirc, offresrabi, stdpn, exczon},

(*measurement and initialisation sequence*)
miseq[q__]:=(Length@{q}>1)&&((Sort[{q}]===Range[0,Max@q])||(Sort[{q}]===Range[Min@q,-1+qubitsnum]));
(* final init state is 011..110 *)
initf[q__]:=Flatten@{
		(* perfect init *)
		Table[
		If[(i===0)||(i===-1+qubitsnum),{Subscript[Damp, i][1.]},{Subscript[Damp, i][1.],Subscript[X, i]}],{i,{q}}],
		(* the errors *)
		Table[
		(* the ends has more measurement+flip errors *)
		If[(pair==={0,1})||(pair==={qubitsnum-2,qubitsnum-1}),
			{Subscript[Depol, Sequence@@pair][erread["01"][[1]]],Subscript[Deph, Sequence@@pair][erread["01"][[1]]]}
		,
		(* the middle has rot errors *)
			{Subscript[Depol, Sequence@@pair][er1xy[Min@pair][[1]]],Subscript[Deph, Sequence@@pair][er1xy[Min@pair][[2]]]}
		],{pair, grouptwo[{q}]}]
		};
		
measf[q__]:=Flatten@{
		(* the errors *)
		Table[
		(* the ends has more measurement+flip errors *)
		If[(pair==={0,1})||(pair==={qubitsnum-2,qubitsnum-1}),
			{Subscript[Depol, Sequence@@pair][erread["01"][[1]]],Subscript[Deph, Sequence@@pair][erread["01"][[1]]]}
		,
		(* the middle has rot errors *)
			{Subscript[Depol, Sequence@@pair][er1xy[Min@pair][[1]]],Subscript[Deph, Sequence@@pair][er1xy[Min@pair][[2]]]}
		],{pair, grouptwo[{q}]}]
		,
		Table[Subscript[M, i],{i,{q}}]
		,
		(* more errors after measurement, non-ideal projected state *)
		Table[
		(* the ends has more measurement+flip errors *)
		If[(pair==={0,1})||(pair==={qubitsnum-2,qubitsnum-1}),
			{Subscript[Depol, Sequence@@pair][erread["01"][[1]]],Subscript[Deph, Sequence@@pair][erread["01"][[1]]]}
		,
		(* the middle has rot errors *)
			{Subscript[Depol, Sequence@@pair][er1xy[Min@pair][[1]]],Subscript[Deph, Sequence@@pair][er1xy[Min@pair][[2]]]}
		],{pair, grouptwo[{q}]}]
		};


stdpn[q_,dd:True]:=If[stdpassivenoise,{Subscript[Deph, q][.5(1-E^(-\[CapitalDelta]t/t2[q]))],Subscript[Depol, q][.75(1-E^(-\[CapitalDelta]t/t1[q]))]},{}];
stdpn[q_,dd:False]:=If[stdpassivenoise,{Subscript[Deph, q][.5(1-E^(-\[CapitalDelta]t/t2s[q]))],Subscript[Depol, q][.75(1-E^(-\[CapitalDelta]t/t1[q]))]},{}];
passivenoisecirc[q_]:=Flatten@{If[q<qubitsnum-1 && AssociationQ@exchangerotoff,Subscript[C, q][Subscript[Rz, q+1][(\[CapitalDelta]t/\[Pi])*exchangerotoff[q]]],{}],stdpn[q,stdpassivenoise]};

(*single rotation noise *)
offresrabi[q_,\[Theta]_]:=If[offresonantrabi,Table[Subscript[U, j][OffResRabiOsc[rabifreq[q],rabifreq[j],Abs[\[Theta]]]],{j,Delete[qubits,q+1]}],{}];

(*Exchange rotation C-Rz[j] interaction when CZ gate on*)
exczon[targ_]:=If[ListQ@exchangeroton,Subscript[C, #-1][Subscript[Rz, #][exchangeroton[[targ,#]]]]&/@Delete[Range[qubitsnum-1],targ],{}];




<|
(*no hidden qubits/ancilla here *)
DeviceDescription -> "Delft Silicon device with "<>ToString[qubitsnum]<>"-qubits arranged as a linear array with nearest-neighbor connectivity and control qubits are the lower ones.",
NumAccessibleQubits -> qubitsnum,
NumTotalQubits -> qubitsnum,

Aliases -> {
	Subscript[Wait, q__][t_] :> {}
	,
	Subscript[Init, i__]:> {}
	}
	,	
Gates ->{
	Subscript[Wait, q__][t_]:><|
		NoisyForm->Flatten@Table[passivenoisecirc[i],{i,Flatten@{q}}],
		GateDuration->t
	|>,
(* Measurements and initialisation *)
	Subscript[M, q__]/; miseq[q] :><|
		NoisyForm-> measf[q],
		GateDuration->durmeas[StringRiffle[Sort@{q},""]]
	|>,
	Subscript[Init, q__]/; miseq[q] :><|
		NoisyForm-> initf[q],
		GateDuration->durinit[StringRiffle[Sort@{q},""]]
	|>,
		
	
(* Singles *)
	Subscript[Rz, q_][\[Theta]_]:> <|
		NoisyForm->{Subscript[Rz, q][\[Theta]]},
		GateDuration->0(*virtual gate, instant, noiseless*)
		|>,
	Subscript[Rx,q_][\[Theta]_]:><|
		NoisyForm->{Subscript[Rx, q][\[Theta]],Subscript[Depol, q][er1xy[q][[1]]],Subscript[Deph, q][er1xy[q][[2]]], Sequence@@offresrabi[q,\[Theta]]},
		GateDuration->Abs[\[Theta]]/(2\[Pi] freqsinglexy[q]) 
	|>,
	Subscript[Ry,q_][\[Theta]_]:><|
		NoisyForm->{Subscript[Ry, q][\[Theta]],Subscript[Depol, q][er1xy[q][[1]]],Subscript[Deph, q][er1xy[q][[2]]], Sequence@@offresrabi[q,\[Theta]]},
		GateDuration->Abs[\[Theta]]/Subscript[(2\[Pi] freqsinglexy[q]), \[Placeholder]]
	|>,
(* Twos *)
	Subscript[C, p_][Subscript[Z, q_]]/; (q-p)===1  :><|
		(*The last bit undo the exchange in the passive noise *)
		NoisyForm->{Subscript[C, p][Subscript[Z, q]],Subscript[Depol, p,q][ercz[p][[1]]],Subscript[Deph, p,q][ercz[p][[2]]],Sequence@@exczon[q]}, 
		GateDuration->0.5/freqcz[p] 
	|>,
	
	(* This is technically single x conditioned on outcome *)
		Subscript[C, p_][Subscript[Rx, q_][\[Theta]_]]/; Abs[(q-p)]===1  :><|
		NoisyForm->{Subscript[C, p][Subscript[Rx, q][\[Theta]]],Sequence@@passivenoisecirc[p],Subscript[Depol, q][er1xy[q][[1]]],Subscript[Deph, q][er1xy[q][[2]]],Sequence@@offresrabi[q,\[Theta]]}, 
		GateDuration->Abs[\[Theta]]/(2\[Pi] freqsinglexy[q])   
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
TrappedIonOxford[OptionsPattern[]]:=With[
{

},

Module[
{\[CapitalDelta]t, miseq,initf,measf, passivenoise},

<|
(*no hidden qubits/ancilla here *)
DeviceDescription -> "Trapped ion device.",
NumAccessibleQubits -> qubitsnum,
NumTotalQubits -> qubitsnum,

Aliases -> {
	Subscript[Wait, q__][t_] :> {}
	,
	Subscript[Init, i__]:> {}
	}
	,	
Gates ->{
	Subscript[Wait, q__][t_]:><|
		NoisyForm->Flatten@Table[passivenoisecirc[i],{i,Flatten@{q}}],
		GateDuration->t
	|>		
},
(* Declare that \[CapitalDelta]t will refer to the duration of the current gate/channel. *)
DurationSymbol -> \[CapitalDelta]t, 
(* Passive noise *)
Qubits :> {
		q_ :> <|
		PassiveNoise ->If[passivenoise===True, passivenoisecirc[q],{}]
		|>		 
		}
		
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
