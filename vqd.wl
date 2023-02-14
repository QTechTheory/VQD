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
SiliconDelft2::usage="Returns device specification of a Silicon device with twice error severity of SiliconDelft.";
SiliconHub::usage="Returns devices specification of a Silicon device based on the device built by the QCSHub.";

(*superconducting qubit devices*)
SuperconductingFZJ::usage="Returns device specification of a Superconducting qubit device based on the device built by Forschungzentrum Juelich.";
SuperconductingHub::usage="Returns device specification of a Superconducting qubit device based on the device built by the QCSHub.";

(*trapped ion devices*)
TrappedIonOxford::usage="Returns device specification of a multi-nodes Trapped ions based on the device built by the Oxford/Hub.";
TrappedIonInnsbruck::usage="Returns device specification of a string of Trapped ions base on the device built by the University of Innsbruck.";

(*rydberg quantum devices/neutral atoms.*)
RydbergHub::usage="Returns device specification of a Rydberg/Neutral Atom device based on the device built by the QCSHub.";
(*RydbergHarvard::usage="Returns device specification of a Rydberg/Neutral Atom device based on the device built in Harvard.";*)
RydbergWisconsin::usage="Returns device specification of a Rydberg/Neutral Atom device based on the device built by the University of Wisconsin.";

(*nuclear-vacancy center devices.*)
NVCenterDelft::usage="Returns device specification of a Nitrogen-Vacancy diamond center device based on the device built by the University of Delft.";
NVCenterHub::usage="Returns device specification of a Nitrogen-Vacancy diamond center device based on the device built by the QCSHub.";

(* toy device *)
ToyDevice::usage="Return a specification with simple standard model.";

ParameterDevices::usage="Show all parameters used to specify all devices. To see parameters related to a device, e.g., NVCenterHub, use Options[NVCenterHub].";

(* General functions  *)
CalcFidelityDensityMatrices::usage="CalcFidelityDensityMatrices[\[Rho],\[Sigma]] fidelity of two density matrices, \[Rho] and \[Sigma] can be density matrix of Quregs. Fidelity of two density matrices.";
PartialTrace::usage="PartialTrace[qureg/density matrix, tracedoutqubits_List]. Return the partial trace as a matrix.";
RandomMixState::usage="RandomMixState[nqubits, nsamples:None]. Return a random mixed quantum density state matrix.";
GenerateOptionTable::usage="GenerateOptionTable[options,columnwidth:{3cm,3cm,10cm}]. Create summary of options in Latex. Very buggy.";
(* Custom gates *)
SWAPLoc::usage="Swap the spatial locations of two qubits";
ShiftLoc::usage="ShiftLoc[v] the physical coordinate of a qubit by some vector v.";
Wait::usage="Wait gate, doing nothing";
Init::usage="Initialise qubit to state |0>";
CZ::usage="Controlled-Z operation";
CRx::usage="Conditional Rx[\[Theta]] rotation on the nuclear 13C NV-center qubit, conditioned on the electron spin state.";
CRy::usage="Conditional Ry[\[Theta]] rotation on the nuclear 13C NV-center qubit, conditioned on the electron spin state.";
Ent::usage="Remote entanglement operation";
Splz::usage="Splz[node, zone_destination]. Split a string of ions in a zone of a trapped-ion Oxford device";
Shutl::usage="Shutl[node,zone_dest]. Shuttle the qubit(s) to the destination zone";
Comb::usage="Comb[node, zone_destination]. Combine a string of ions to a zone of a trapped-ion Oxford device";
PSW::usage="PSW[\[Theta]], parameterised swaps";
UG::usage="Single unitary gate rotation obtained by the evolution of driven qubit, e.g., by laser: UG[\[Phi],\[CapitalDelta],t,\[CapitalOmega]], where \[Phi] is laser phase, \[CapitalDelta] is detuning, t is laser duration, and \[CapitalOmega] is the Rabi frequency.";
SRot::usage="Single qubit gate in a driven Rydberg qubit via two-photon Raman transition: SRot[\[Phi],\[CapitalDelta],t] where \[Phi] is laser phase, \[CapitalDelta] is detuning, t is laser duration.";

(*Visualisations*)
DeviceType::usage="The type of device. Normally, the name of the function that generates it.";
DrawIons::usage="Draw the current string of ions";

(** PlotAtoms options for Rydberg device **)
PlotAtoms::usage="PlotAtoms[rydberg_device]. Plot the atoms of a Rydberg device. Set ShowBlockade->{qubits} to show the blockade radius of a set of qubits. Set ShowLossAtoms->True, to show the atoms that are loss as well. It will be indicated with grey color.
Plot atoms also receive the options of Graphic function.";
PlotAtoms::error="`1`";
Options[PlotAtoms]={
ShowBlockade->{},
ShowLossAtoms->False
};
ShowBlockade::usage="List the qubits to draw the blockade radius.";
ShowLossAtoms::usage="Set true to show the atoms lost into the evironment. This shows the last coordinate before being lost.";
(* Options for functions *)
Parallel::usage="Parallel options in arrangement of gates: False, Default, All. 
False: serial, default: parallel according to the device specification, and All: full quantum parallel";
MapQubits::usage="Options in the CircTrappedIons[] that maps the local qubits into the total qubits of the total (large) density matrix.";

(*Circuits arrangement according to devices*)
Options[CircTrappedIons]={MapQubits->True, Parallel->False};
Options[CircSiliconDelft]={Parallel->False};
Options[CircRydbergHub]={Parallel->False};

CircTrappedIons::usage="CircTrappedIons[circuit, device, MapQubits->True, Parallel->False]. Circuit arrangement according to the device. Note that Parallle->True is not available yet.";
CircSiliconDelft::usage="CircSiliconDelft[circuit, device, Parallel->False]. Circuit arrangement according to the device. Note that Parallle->True is not available yet.";
CircRydbergHub::usage="CircRydbergHub[circuit, device, Parallel->(False, True)]";
Serialize::usage="Serialize circuit. Every quantum operation is done without concurency.";
CircTrappedIons::error="`1`";
CircSiliconDelft::error="`1`";
CircRydbergHub::error="`1`";
Serialize::error="`1`";
BeginPackage["`ParameterDevices`"];
(* Parameters *)
AtomLocations::usage="Three-dimensional physical locations of each atom/qubit.";
Anharmonicity::usage="The anharmonicity in the Superconducting device capturing the capacitor property.";
BFProb::usage="Probability of bit-flip error";
BField::usage="The electromagnetic field strength in the z-direction from the lab reference with unit Tesla.";
BlockadeRadius::usage="Short-range dipole-dipole interaction of Rydberg atoms in \[Mu]s. This allows multi-qubit gates.";
Connectivity::usage="Show the connectivity graph of a Superconducting qubit device, where the arrow show possible direction of the cross-resonant ZX gates.";
DurTwoGate::usage="Duration of two qubit gates with rotation of \[Pi]";
DurMeas::usage="Duration of measurement";
DurInit::usage="Duration of initialisation";
DurRead::usage="Readout duration in \[Mu]s";
DurShuffle::usage="Duration to shuffle location of ions";
DurMove::usage="Duration for physically moving operation in Trapped Ions such as Splz and Comb.";
DurRxRy::usage="Duration to run the rotation gates Rx and Ry, that is fixed regardless the angle.";
DurZX::usage="Duration of the resonance ZX gate on the superconducting qubits that is fixed regardless the angle.";
DurZZ::usage="Duration of the siZZle ZZ gate on the superconducting qubits that is fixed regardless the angle.";
DDActive::usage="Apply dynamical decoupling: use T2 in the model if set True, otherwise use T2* if set False.";
entangling::usage="Crosstalk error model pre equation (4) on applying the XX M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
ErrCT::usage="Error coefficient of the crosstalk with entanglement model";
ErrSS::usage="Error coefficient of the crosstalk with stark shift model";
EFSingleXY::usage="Error fraction/ratio, {depolarising, dephasing} of the single qubit X and Y rotations. Sum of the ratio must be 1 or 0 (off).";
EFSingle::usage="Error fraction/ratio, {depolarising, dephasing} of the single qubit X, Y, and Z rotations. Sum of the ratio must be 1 or 0 (off).";
EFTwo::usage="Error fraction/ratio, {depolarising, dephasing} of controlled-rotation. Sum of the ratio must be 1 or 0 (off).";
EFCZ::usage="Error fraction/ratio, {depolarising, dephasing} of the controlled-Z gates. Sum of the ratio must be 1 or 0 (off).";
EFInit::usage="Error fraction/ratio {depolarising, dephasing} of the initialisation gate.  Sum of the ratio must be 1 or 0 (off).";
EFCRot::usage="Error fraction/ratio, {depolarising, dephasing} of the controlled-Rx and -Ry gates. Sum of the ratio must be 1 or 0 (off).";
EFEnt::usage="Error fraction/ratio {depolarising, dephasing} of remote entanglement.  Sum of the ratio must be 1 or 0 (off).";
EFRead::usage="Error fraction/ratio of {depolarising,dephasing} of the readout. Sum of the ratio must be 1 or 0 (off)." ;
ExchangeRotOn::usage="Maximum interaction j on the passive qubit crosstalk when applying CZ gates; The noise form is C[Rz[j.\[Theta]]] It must be a square matrix with size (nqubit-2)x(nqubit-2).";
ExchangeRotOff::usage="Crosstalks error C-Rz[ex] on the passive qubits when not applying two-qubit gates.";
ExcitedInit::usage="The probability/fraction of the population excited in the thermal state. This is the same initialisation state.";
ExchangeCoupling::usage="The exchange coupling strength of resonators in Superconducting devices";
FidCRot::usage="Fidelity of conditional rotation in NV-center obtained by dynamical decoupling and RF pulse.";
FidSingleXY::usage="Fidelity(ies) of single Rx[\[Theta]] and Ry[\[Theta]] rotations obtained by random benchmarking.";
FidSingleZ::usage="Fidelity(ies) of single Rz[\[Theta]] rotation obtained by random benchmarking.";
FidSingle::usage="Fidelity(ies) of single rotations: Rx[\[Theta]], Ry[\[Theta]], Rz[\[Theta]] obtained by random benchmarking.";
FidTwo::usage="Fidelity(ies) of two qubit gates obtained by random benchmarking.";
FidEnt::usage="Fidelity of remote entanglement operation.";
FidMeas::usage="Fidelity of measurement";
FidInit::usage="Fidelity of qubit initialisation";
FidCZ::usage="Fidelity(ies) of the CZ gates.";
FreqSingleXY::usage="Rabi frequency(ies) for the single X- and Y- rotations with unit MHz";
FreqSingleZ::usage="Rabi frequency(ies) for the single Z- rotations with unit MHz";
FreqCZ::usage="Rabi frequency(ies) for the CZ gate with unit MHz.";
FreqEnt::usage="Frequency of remote entanglement.";
FreqCRot::usage="Frequency of conditional rotation in NV-center obtained by dynamical decoupling and RF pulse.";
FreqWeakZZ::usage="Frequency of coherent cross-talk noise in form of ZZ-coupling that slowly entangle the qubits.";
FidRead::usage="Readout fidelity";
GlobalField::usage="Global magnetic field fluctuation in NV-center due to C13 bath. It causes dephasing on the electron in 4\[Mu]s.";
GlobalFieldDetuning::usage="TODO:? forgot";
LossAtoms::usage="Device key in the RydbergHub device that identifies atoms lost to the environment.";
LossAtomsProbability::usage="Device key in the RydbergHub that identifies probability of atoms lost to the environment due to repeated measurement.";
MSCrossTalk::usage="(entangling OR starkshift) The crosstalk model in applying M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
MoveSpeed::usage="The speed of moving atom \[Mu]m/\[Mu]s in Neutral Atom systems.";
Nodes::usage="Entire nodes of a trapped ions system <|node1 -> number_of_qubits_1, ... |>";
NIons::usage="The total number of ions in a trapped ion device.";
OffResonantRabi::usage="Put the noise due to off-resonant Rabi oscillation when applying single qubit rotations.";
qubitsNum::usage="The number of physical active qubits for computations.";
QubitFreq::usage="The Qubit frequency for each qubit with unit MHz.";
QMap::usage="Show maps from nodes in trapped ions to the actual emulated qubits";
Meas::usage="Perform measurement on the qubits";
ProbLeakInit::usage="Leakage probability in the Rydberg initialisation. The noise is decribed with non-trace-preserving map.";
ProbLeakCZ::usage="Leakage probability in executi multi-controlled-Z.";
ProbLossMeas::usage="Probability of phyiscal atom loss due to measurement.";
ProbLeakMove::usage="Probability of leakage in the process of moving atoms.";
ProbBFRot::usage="Assymetric Bit-flip probability on single rotation operation. {01->p1, 10->p2}";
RydbergRabiFreq::usage="Rydberg Rabi frequency";
RabiFreq::usage="The Rabi frequency frequency in average or on each qubit with unit MHz.";
RepeatRead::usage="The number of repeated readout (n) performed. The final fidelity of readout is \!\(\*SuperscriptBox[\(FidRead\), \(n\)]\).";
starkshift::usage="Crosstalk error model using stark shift on modeled with MS gate Exp[-i\[Theta]XX], M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
ScatProb::usage="Scattering probability";
ShowNodes::usage="Draw all Ions on every nodes within the zones";
StdPassiveNoise::usage="Set to True/False. Use the standard passive noise that involves T1, T2 or T2s inputs.";
T1::usage="T1 duration(s) in \[Mu]s. Exponential decay time for the state to be complete mixed.";
T2::usage="T2 duration(s) in \[Mu]s. Exponential decay time for the state to be classical with echo applied.";
T2s::usage="T2* duration(s) in \[Mu]s. Exponential decay time for the state to be classical.";
TwoGateFreq::usage="Resonant frequency of two qubit gates";
UnitLattice::usage="The unit lattice AtomLocations in \[Mu]s. This gives access to internal device parameter in RydbergHub device.";
VacLifeTime::usage="The lifetime of the qubit array is limited by its vacuum lifetime, where T1=VacTime/Nqubits.";
ZZ::usage="The siZZle gate on a Superconducting device. Implemented by Exp[-(i\[Theta]/2) ZZ].";
ZX::usage="The cross-resonance gate on a Superconducting device. Implemented by Exp[-(i\[Theta]/2) ZX].";
ZZPassiveNoise::usage="The swtich for ZZ interaction passive noise in the Superconducting device.";
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
FidSingle::error="`1`";
FidTwo::error="`1`";
FidMeas::error="`1`";
FidInit::error="`1`";
FreqCZ::error="`1`";
QubitFreq::error="`1`";
Nodes::error="`1`";
qubitsNum::error="`1`";
OffResonantRabi::error="`1`";
ProbLeakInit::error="`1`";
ProbLeakCZ::error="`1`";
ProbLossMeas::error="`1`";
ProbLeakMove::error="`1`";
RabiFreq::error="`1`";
StdPassiveNoise::error="`1`";
T1::error="`1`";
T2::error="`1`";
T2s::error="`1`";
FidCRotXY::warning="`1`";
FidCZ::warning="`1`";
FidSingleXY::warning="`1`";
FidSingle::warning="`1`";
FidTwo::warning="`1`";
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

Subscript[UG, q_][\[Phi]_,\[CapitalDelta]_,t_,\[CapitalOmega]_]:=Module[{v\[CapitalDelta],vt,\[Omega]},
	Subscript[U, q][FullSimplify[
	{{Cos[\[Omega] vt/2]-I v\[CapitalDelta]/\[Omega] Sin[\[Omega] vt/2],-I \[CapitalOmega]/\[Omega] Sin[\[Omega]  vt/2] E^(I \[Phi])},
	 {-I \[CapitalOmega]/\[Omega] Sin[\[Omega] vt/2]E^(I \[Phi]),     Cos[\[Omega] vt/2]+I \[CapitalDelta]/\[Omega] Sin[\[Omega] vt/2]}
	}//.{v\[CapitalDelta]->\[CapitalDelta],\[Omega]->Sqrt[\[CapitalOmega]^2+v\[CapitalDelta]^2],vt->t},{\[CapitalOmega]>=0,v\[CapitalDelta]>=0}]
	]
]

fid2DepolDeph::usage = "fid2DepolDeph[totFid, {ratio.Depol,ratio.Deph}, nqerr, err, avgfid:True]. 
Return the parameter for depolarizing and dephasing noise that gives the total fidelity totFid [0,1], where totFid is the average fidelity obtained from random benchmarking (set avgfid:False if it's the worst fidelity)
correction is a constant adjusted to params of dephasing to get average fidelity from the worst fidelity.
";
fid2DepolDeph[totfid_, errratio_, nqerr_, errval_, avgfid_:True] := Module[
  {sol, rate, pdepol, pdeph, entfid,d,p},
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
entfid2DepolDeph[entfid_,errratio_,errval_:FidEnt]:=Module[{pdepol,pdeph,sol,rate},
(** Obtain the error parameters given entanglement fidelity. This works only for the trapped ion network fidelity.
	This is obtained by calculating <\[CapitalPsi]^+|Subscript[Depol, 0,1](d)Subscript[Deph, 0,1](p)|\!\(
\*SuperscriptBox[\(\[CapitalPsi]\), \(+\)]
\*SubscriptBox[\(>\), \(0, 1\)]\)
**)
	sol=Solve[1-rate*errratio[[1]]*4/5-rate*errratio[[2]]*2/3+rate^2*errratio[[1]]*errratio[[2]]*32/45==entfid,{rate}];
	{pdepol, pdeph}= errratio*Min@Abs[rate/.sol];
	If[pdepol>15/16,
		Message[errval::warning,StringForm["(warning) Fidelity might be too low; 2-qubit depolarization parameter is ``. Set it to 15/16.",pdepol]];
		pdepol=15/16;
	];
	If[pdeph>3/4,
		Message[errval::warning,StringForm["(warning) Fidelity might be too low; 2-qubit dephasing parameter is ``. Set it to 3/4.",pdeph]];
		pdeph=3/4;
	];
	{pdepol,pdeph}
]
grouptwo::usage="grouptwo[list], group a list into two elements";
grouptwo[list_]:=ReplaceList[Sort@list,{p___,a_,b_,q___}:>{a,b}]
(**** EXTRA_QUANTUM_CHANNELS ****)
bitFlip1::usage="bitFlip1[] Return the list of kraus operators with 1- or 2- qubit bitflip error";
bitFlip1[fid_]:=With[{e=1-fid},
	{Sqrt[1-e]*{{1,0},{0,1}},Sqrt[e]*{{0,1},{1,0}}}
	]
bitFlip2[fid_]:=With[{e=1-fid},{
		Sqrt[1-e]*IdentityMatrix[4],
		Sqrt[e/3]*{{0,1,0,0},{1,0,0,0},{0,0,0,1},{0,0,1,0}},
		Sqrt[e/3]*{{0,0,1,0},{0,0,0,1},{1,0,0,0},{0,1,0,0}},
		Sqrt[e/3]*{{0,0,0,1},{0,0,1,0},{0,1,0,0},{1,0,0,0}}}]
(* 
Generalised amplitude damping (Nielsen&Chuang, P.382).
Common description of T1 decay, by setting \[Gamma](t)=1-\[ExponentialE]^(-t/T1),
 which describes the shrinking Bloch sphere into the ground state |0\[RightAngleBracket],
with probability p.
*)
(* Subscript[\[Rho], \[Infinity]]=p|0X0|+(1-p)|1X1| *)
Subscript[gAmp, q_][\[Gamma]_,p_]:=Subscript[Kraus, q][{
	Sqrt[p]*{{1,0},{0,Sqrt[1-\[Gamma]]}},
	Sqrt[p]*{{0,Sqrt[\[Gamma]]},{0,0}},
	Sequence@@If[p<1, 
	{Sqrt[1-p]*{{Sqrt[1-\[Gamma]],0},{0,1}},
	Sqrt[1-p]*{{0,0},{Sqrt[\[Gamma]],0}}},{}]}]
				
								
(***** TOY_DEVICE *****)
ToyDevice[OptionsPattern[]]:=With[
	{
	qubitsnum=OptionValue[qubitsNum],
	t1=OptionValue[T1],
	t2=OptionValue[T2],
	stdpassivenoise=OptionValue[StdPassiveNoise],
	fidsingle=OptionValue[FidSingle],
	fidtwo=OptionValue[FidTwo],
	efsingle=OptionValue[EFSingle],
	eftwo=OptionValue[EFTwo],
	rabifreq=OptionValue[RabiFreq],
	twogatefreq=OptionValue[TwoGateFreq]
	},
	
	Module[
	{\[CapitalDelta]T,erone,ertwo},
	erone=fid2DepolDeph[fidsingle,efsingle,1,FidSingle,True];
	ertwo=fid2DepolDeph[fidtwo,eftwo,2,FidTwo,True];
	<|
	(*no hidden qubits/ancilla here *)
	DeviceType->"Toy",
	DeviceDescription -> "Toy device with "<>ToString[qubitsnum]<>"-qubits arranged as a linear array with nearest-neighbor connectivity.",
	NumAccessibleQubits -> qubitsnum,
	NumTotalQubits -> qubitsnum,
	Aliases -> {
	Subscript[PSW, p_,q_][\[Theta]_]:>Subscript[U, p,q][{{1,0,0,0},{0,E^((I \[Theta])/2) Cos[\[Theta]/2],-I E^((I \[Theta])/2) Sin[\[Theta]/2],0 },{0,-I E^((I \[Theta])/2) Sin[\[Theta]/2], E^((I \[Theta])/2) Cos[\[Theta]/2],0 },{0,0,0,1}}]  
	},	
	Gates ->{	
	(* Singles *)
		Subscript[Rx,q_][\[Theta]_]:><|
			NoisyForm->Flatten@{Subscript[Rx, q][\[Theta]],Subscript[Depol, q][erone[[1]]],Subscript[Deph, q][erone[[2]]]},
			GateDuration->Abs[\[Theta]]/rabifreq
		|>,
		Subscript[Ry,q_][\[Theta]_]:><|
			NoisyForm->Flatten@{Subscript[Ry, q][\[Theta]],Subscript[Depol, q][erone[[1]]],Subscript[Deph, q][erone[[2]]]},
			GateDuration->Abs[\[Theta]]/rabifreq
		|>,
			Subscript[Rz,q_][\[Theta]_]:><|
			NoisyForm->Flatten@{Subscript[Rz, q][\[Theta]],Subscript[Depol, q][erone[[1]]],Subscript[Deph, q][erone[[2]]]},
			GateDuration->Abs[\[Theta]]/rabifreq
		|>,
	(* Twos *)
			Subscript[C, p_][Subscript[Rx, q_][\[Theta]_]]/; Abs[q-p]===1  :><|
			NoisyForm->{Subscript[C, p][Subscript[Rx, q][\[Theta]]],Subscript[Depol, p,q][Min[ertwo[[1]]*Abs[\[Theta]/\[Pi]],15/16]],Subscript[Deph, p,q][Min[ertwo[[2]]*Abs[\[Theta]/\[Pi]],3/4]]},
			GateDuration->Abs[\[Theta]]/twogatefreq
		|>,
			Subscript[C, p_][Subscript[Ry, q_][\[Theta]_]]/; Abs[q-p]===1  :><|
			NoisyForm->{Subscript[C, p][Subscript[Ry, q][\[Theta]]],Subscript[Depol, p,q][Min[ertwo[[1]]*Abs[\[Theta]/\[Pi]],15/16]],Subscript[Deph, p,q][Min[ertwo[[2]]*Abs[\[Theta]/\[Pi]],3/4]]},
			GateDuration->Abs[\[Theta]]/twogatefreq
		|>,
			Subscript[C, p_][Subscript[Rz, q_][\[Theta]_]]/; Abs[q-p]===1  :><|
			NoisyForm->{Subscript[C, p][Subscript[Rz, q][\[Theta]]],Subscript[Depol, p,q][Min[ertwo[[1]]*Abs[\[Theta]/\[Pi]],15/16]],Subscript[Deph, p,q][Min[ertwo[[2]]*Abs[\[Theta]/\[Pi]],3/4]]},
			GateDuration->Abs[\[Theta]]/twogatefreq
		|>	
		,(* parameterised swap*)
			Subscript[PSW, p_,q_][\[Theta]_]/; Abs[q-p]===1  :><|
			NoisyForm->{Subscript[PSW, p,q][\[Theta]],Subscript[Depol, p,q][ertwo[[1]]],Subscript[Deph, p,q][ertwo[[2]]]},
			GateDuration->Abs[\[Theta]]/twogatefreq
		|>	
	},
	(* Declare that \[CapitalDelta]t will refer to the duration of the current gate/channel. *)
	DurationSymbol -> \[CapitalDelta]T, 
	(* Passive noise *)
	Qubits :> {
			q_Integer :> <|
			PassiveNoise ->If[stdpassivenoise,{Subscript[Depol, q][0.75(1-E^(-\[CapitalDelta]T/t1))],Subscript[Deph, q][0.5(1-E^(-\[CapitalDelta]T/t2))]},{}]	
			|>		 
			}	
		|>
	]
]
(******************************* SUPERCONDUCTING_QUBITS *********************************)
(*Returns graph that show connectivity*)
graphConnectivity[nqubit_,coupling_,freq_]:=Module[{g,dedge},
	(*directed edges, control->target*)
	dedge=If[freq[#[[1]]]>freq[#[[2]]],#[[1]]->#[[2]],#[[2]]->#[[1]]]&/@Keys[coupling];
	(* Basic graph *)
	g=Graph[Range[0,nqubit-1],dedge,VertexWeight->Values@freq,EdgeWeight->Values@coupling];
	DirectedGraph[g,VertexSize->0.5,BaseStyle->{15,Bold,FontFamily->"Serif"},
	VertexLabels->{v_:>Placed[{AnnotationValue[{g,v},VertexWeight],"Q"<>ToString[v]},{Center,{After,Above}}]},
	EdgeLabels->"EdgeWeight",
	GraphLayout->"SpringEmbedding",EdgeStyle->{Thick},VertexStyle->{Yellow,EdgeForm[None]}]
]
SuperconductingHub[OptionsPattern[]]:=With[
	{
	qubitsnum=OptionValue[qubitsNum],
	t1=OptionValue[T1],
	t2=OptionValue[T2],
	excitedinit=OptionValue[ExcitedInit],
	qubitfreq=OptionValue[QubitFreq],
	exchangecoupling=OptionValue[ExchangeCoupling],
	anharmonicity=OptionValue[Anharmonicity],
	fidread=OptionValue[FidRead],
	durmeas=OptionValue[DurMeas],
	durrxry=OptionValue[DurRxRy],
	durzx=OptionValue[DurZX],
	durzz=OptionValue[DurZZ],
	stdpassivenoise=OptionValue[StdPassiveNoise],
	zzpassivenoise=OptionValue[ZZPassiveNoise]
	},
	(* Assertions *)
	Catch@If[CountDistinct[Values@qubitfreq]==Length@qubitfreq, Throw@Message[QubitFreq::error,"All qubits frequencies must be distinct. Fix the QubitFreq value."]];
	Module[
	{ccv,ug, stdpn, zz, zzpn, lessNeighbor, zzon, passivenoise, \[CapitalDelta]t, activeq,counter=0,init=True},
	ccv=graphConnectivity[qubitsnum,exchangecoupling,qubitfreq];	
	(* undirected graph *)
	ug=UndirectedGraph[ccv];
	activeq=<|Table[q->False,{q,Range[0,qubitsnum-1]}]|>;
	(* if ZZ gate is on *)
	zzon=False;
	(* Free-induction T1- T2- passive noise decays *)
	Subscript[stdpn, q_][\[CapitalDelta]t_]:=If[stdpassivenoise&&\[Not]activeq[q]&&\[CapitalDelta]t>0,
					{Subscript[gAmp, q][(1-E^(-\[CapitalDelta]t/t1[q])),1-excitedinit[q]],
					Subscript[Deph, q][0.5(1-E^(-\[CapitalDelta]t/t2[q]))]},
					{}];
	(* Fixed ZZ-interaction on passive noise, where p is the first/control and q is the later/target *)
	Subscript[zz, c_, t_]:=With[
	{\[CapitalDelta]ct=Abs[qubitfreq[c]-qubitfreq[t]],\[Alpha]t=anharmonicity[t],\[Alpha]c=anharmonicity[c],J=If[KeyExistsQ[exchangecoupling,c\[UndirectedEdge]t],exchangecoupling[c\[UndirectedEdge]t],exchangecoupling[t\[UndirectedEdge]c]]},
			R[J^2 (1/(\[CapitalDelta]ct-\[Alpha]t)-1/(\[CapitalDelta]ct+\[Alpha]t)),Subscript[Z, c] Subscript[Z, t]]];
	(* list neighbor qubits with ordering less than q in graph g. This module is used in the superconducting device. *)
	lessNeighbor[q_]:=List@@@Select[EdgeList[NeighborhoodGraph[ccv,q]],q>First@DeleteElements[List@@#,{q}]&];
	Subscript[zzpn, q_][\[CapitalDelta]t_]:=Module[{ng=lessNeighbor[q],noise},
			noise=If[zzpassivenoise, 
			If[\[Not]zzon,(If[\[CapitalDelta]t>0&&And[\[Not]activeq[#[[1]]],\[Not]activeq[#[[2]]]],Subscript[zz, #[[1]],#[[2]]],{}])&/@ng,{}],{}];
			noise
			];
			
			
	<|
	(*no hidden qubits/ancilla here *)
	DeviceType->"Superconducting",
	DeviceDescription -> ToString[qubitsnum]<>"-qubit of Superconducting transmon qubits based on Josephson junctions",
	NumAccessibleQubits -> qubitsnum,
	NumTotalQubits -> qubitsnum,
	Connectivity->ccv,
	Aliases -> {
		Subscript[ZZ, p_,q_]:>R[\[Pi],Subscript[Z, p] Subscript[Z, q]],
		Subscript[ZX, p_,q_]:>R[\[Pi],Subscript[Z, p] Subscript[X, q]],
			(* Subscript[\[Rho], \[Infinity]]=p|0X0|+(1-p)|1X1| *)
		Subscript[Init, q_]:>Sequence@@{},
		Subscript[Wait, q_]:>Sequence@@{}
	},	
	Gates ->{
	(* Initialisation and measurement *)
	Subscript[Init, q_]/;(init==True):><|
		NoisyForm->{Subscript[Init, q],Subscript[gAmp, q][1,1-excitedinit[q]]},
		GateDuration->0
	|>,
	(* should be placed at the end only, but I don't know how :-/ *)
	Subscript[Read, q_]:><|
		(* mess up the final result as well *)
		NoisyForm->{Subscript[Depol, q][1-fidread[q]],Subscript[M, q],Subscript[Depol, q][1-fidread[q]]}
	
	|>,
	(* doing nothing, is equivalent to being passive *)
	Subscript[Wait, q_][t_]:><|
		NoisyForm->Flatten@{Subscript[stdpn, q][\[CapitalDelta]t],Subscript[zzpn, q][\[CapitalDelta]t]},
		GateDuration->t,
		UpdateVariables->Function[
			activeq[q]=False
			]
	|>,
	(* Singles: commenting the angle rule for VQE program *)
		Subscript[Rx,q_][\[Theta]_](*/;And[(-\[Pi]<=\[Theta]<=\[Pi]),Abs[\[Theta]]>0]*):><|
			NoisyForm->{Subscript[Rx, q][\[Theta]]},
			GateDuration->durrxry,
			UpdateVariables->Function[
			activeq[q]=True;
			init=False;
			]
		|>,
		Subscript[Ry,q_][\[Theta]_](*/;And[(-\[Pi]<=\[Theta]<=\[Pi]),Abs[\[Theta]]>0]*):><|
			NoisyForm->{Subscript[Ry, q][\[Theta]]},
			GateDuration->durrxry,
			UpdateVariables->Function[
			activeq[q]=True;
			init=False;
			]
		|>,
		
		Subscript[Rz,q_][\[Theta]_]:><|
			NoisyForm->Flatten@{Subscript[Rz, q][\[Theta]]},
			GateDuration->0.
		|>,
		(* siZZle gate *)
		Subscript[ZZ, p_,q_]/;EdgeQ[ug,p\[UndirectedEdge]q] :><|
			NoisyForm->{Subscript[ZZ, p,q]},
			GateDuration->durzz,
			UpdateVariables->Function[
				init=False;
				zzon=True;	
				activeq[q]=True;
				activeq[p]=True;
			]
		|>,
		Subscript[ZX, p_,q_]/;EdgeQ[ccv,p\[DirectedEdge]q] :><|
			NoisyForm->{Subscript[ZX, p,q]},
			GateDuration->durzx,
			UpdateVariables->Function[
				init=False;
				activeq[q]=True;
				activeq[p]=True;
			]
		|>
		
	},
	(* Passive noise *)
	(* Declare that \[CapitalDelta]t will refer to the duration of the current gate/channel. *)
	DurationSymbol -> \[CapitalDelta]t, 
	Qubits :> {
		q_ :> <|
				(* reset the zzon here *)
				PassiveNoise -> Flatten@{Subscript[stdpn, q][\[CapitalDelta]t],Subscript[zzpn, q][\[CapitalDelta]t]}
				,
				UpdateVariables->Function[
				counter++; 
				If[counter==qubitsnum,	
					zzon=False;
					Table[activeq[j]=False,{j,Keys@activeq}];
					counter=0;
				]
				]	
				|>	
				
			}
		|>
	]
]

(**************************** NVCENTER_DELFT *****************************)
NVCenterDelft[OptionsPattern[]]:=With[
{
	qubitsnum=OptionValue@qubitsNum,
	t1=OptionValue@T1,
	t2=OptionValue@T2,
	freqcrot=OptionValue@FreqCRot,
	freqsinglexy=OptionValue@FreqSingleXY,
	freqsinglez=OptionValue@FreqSingleZ,
	fidcrot=OptionValue@FidCRot,
	fidsinglexy=OptionValue@FidSingleXY,
	fidsinglez=OptionValue@FidSingleZ,
	efsinglexy=OptionValue@EFSingleXY,
	efcrot=OptionValue@EFCRot,
	fidinit=OptionValue@FidInit,
	fidmeas=OptionValue@FidMeas,
	durmeas=OptionValue@DurMeas,
	durinit=OptionValue@DurInit,
	freqweakzz=OptionValue@FreqWeakZZ
},


Module[{\[CapitalDelta]t, stdpn,pn,ersinglexy,ersinglez,ercrot,nuclearq,weakzz},
(** standard passive noise **) 
	stdpn[q_,dur_]:={Subscript[Depol, q][0.75(1-E^(-dur/t1[q]))],Subscript[Deph, q][0.5(1-E^(-dur/t2[q]))]};
	(* 
	some note: cross-talk ZZ-coupling in order of Hz on passive noise, Phase entanglement among nuclear spins 
	Subscript[C, n1][Subscript[Rz, n2][\[CapitalDelta]t]], weak rotation for all combinations of n1,n2, few Hz.
	implement: Exp[-i dur ZZ]
	*)	
	weakzz[q_,dur_]:=R[dur*\[Pi]/freqweakzz,Subscript[Z, #] Subscript[Z, #2]]&@@@Subsets[DeleteElements[Range[1,qubitsnum-1],{q}],{2}];

	pn[q_,dur_]:=If[NumberQ@freqweakzz,Flatten@{stdpn[q,dur],weakzz[q,dur]},stdpn[q,dur]];


(* error parameters *)
	ersinglexy=fid2DepolDeph[#,efsinglexy,1,FidSingleXY]&/@fidsinglexy;
	ersinglez=fid2DepolDeph[#,{0,1},1,FidSingleZ]&/@fidsinglez;
	ercrot=fid2DepolDeph[#,efcrot,2,FidCRot]&/@fidcrot;

<|
	(* A helpful description of the device *)
	DeviceDescription ->"One node of an NV center, where qubit 0 is the electronic spin. It has start connectivity with qubit 0 at the center.",
	(* The number of accessible qubits. This informs the qubits that a user's circuit can target *)
	NumAccessibleQubits -> qubitsnum,	
	(* The total number of qubits which would be needed to simulate a circuit prescribed by this device.
	 * This is > NumAccessibleQubits only when the spec uses hidden qubits for advanced noise modelling *)
	NumTotalQubits -> qubitsnum,

	(* Aliases are useful for declaring custom events. At this stage the specification is noise-free (but see later) *)
Aliases -> {
	Subscript[Init, q_]:> Sequence@@{},
	Subscript[Wait, q__][t_]:>Sequence@@{},
	Subscript[CRx, e_,n_][\[Theta]_]:>Sequence@@{Subscript[U, e,n][{{Cos[\[Theta]/2],0,-I Sin[\[Theta]/2],0},{0,Cos[\[Theta]/2],0,I Sin[\[Theta]/2]},{-I Sin[\[Theta]/2],0,Cos[\[Theta]/2],0},{0,I Sin[\[Theta]/2],0,Cos[\[Theta]/2]}}]},
	Subscript[CRy, e_,n_][\[Theta]_]:>Sequence@@{Subscript[U, e,n][{{Cos[\[Theta]/2],0,-Sin[\[Theta]/2],0},{0,Cos[\[Theta]/2],0,Sin[\[Theta]/2]},{Sin[\[Theta]/2],0,Cos[\[Theta]/2],0},{0,-Sin[\[Theta]/2],0,Cos[\[Theta]/2]}}]}
	},	
Gates ->{
(** exclusively on ELECTRON SPIN, q===0 **)
		Subscript[Init, q_]/;q===0 :> <|
			NoisyForm -> Circuit[Subscript[Damp, q][fidinit]], 
			GateDuration -> durinit
		|>,
		Subscript[M, q_]/;q===0 :> <|
			NoisyForm -> {Subscript[X, 0],Subscript[Damp, 0][1-fidmeas],Subscript[X, 0],Subscript[M, 0]},
			GateDuration -> durmeas
		|>,
		(* Electron and nuclear spins *)
		(* A simple 'wait' instruction, useful for padding a circuit. The content inside the table should be the same as the passive noise section below (copy / paste it) *)
		Subscript[Wait, qubits__][dur_]:> <|
			NoisyForm->stdpn[#,dur]&/@Flatten{qubits},
			GateDuration->dur
		|>,
		
		(* Z-rotations are unconditional for both electron and nuclear spins *)
		Subscript[Rz, q_][\[Theta]_] :> <|
			NoisyForm -> Circuit[Subscript[Rz, q][\[Theta]]Subscript[Deph, q][Min[ersinglez[q][[2]]*Abs[\[Theta]]/\[Pi],0.5]]],
			GateDuration -> Abs[\[Theta]]/freqsinglez[q]
		|>, 
		
		Subscript[Rx, q_][\[Theta]_]/;q===0 :> <|
			NoisyForm -> Circuit[Subscript[Rx, q][\[Theta]]Subscript[Depol, q][Min[ersinglexy[q][[1]]*Abs[\[Theta]]/\[Pi],0.75]]Subscript[Deph,q][Min[ersinglexy[q][[2]]Abs[\[Theta]]/\[Pi],0.5]]],
			GateDuration -> Abs[\[Theta]]/freqsinglexy[q]
		|>,
		
		Subscript[Ry, q_][\[Theta]_]/;q===0 :> <|
			NoisyForm -> Circuit[Subscript[Ry, q][\[Theta]]Subscript[Depol, q][Min[ersinglexy[q][[1]]Abs[\[Theta]]/\[Pi],0.75]]Subscript[Deph, q][Min[ersinglexy[q][[2]]Abs[\[Theta]]/\[Pi],0.5]]],
			GateDuration -> Abs[\[Theta]]/freqsinglexy[q]
		|>,
		(* ROTATIONS CONDITIONED ON ELETRON SPIN. electron needs to set at ms=-1: horrible*)
		Subscript[Rx, q_][\[Theta]_]/;q>0 :> <|
			NoisyForm -> Circuit[Subscript[Rx,q][\[Theta]]Subscript[Depol, q][Min[ersinglexy[q][[1]]Abs[\[Theta]]/\[Pi],0.75]]Subscript[Deph, q][Min[ersinglexy[q][[2]]*Abs[\[Theta]]/\[Pi],0.5]]],
			GateDuration -> Abs[\[Theta]]/freqsinglexy[q]
		|>,
		Subscript[Ry, q_][\[Theta]_]/;q>0 :> <|
			NoisyForm -> Circuit[Subscript[Ry,q][\[Theta]]Subscript[Depol, q][Min[ersinglexy[q][[1]]Abs[\[Theta]]/\[Pi],0.75]]Subscript[Deph, q][Min[ersinglexy[q][[2]]*Abs[\[Theta]]/\[Pi],0.5]]],
			GateDuration -> Abs[\[Theta]]/(freqsinglexy[q])
		|>,
		(* Conditional rotations *)
		Subscript[CRx, e_,n_][\[Theta]_] /; (e==0 && n>0):> <|
			(* Its noisy form depolarises the control and target qubits *)
			NoisyForm -> Circuit[Subscript[CRx, e,n][\[Theta]]Subscript[Depol,e,n][Min[ercrot[n][[1]]Abs[\[Theta]]/\[Pi],15/16]]Subscript[Deph,e,n][Min[ercrot[n][[2]]Abs[\[Theta]]/\[Pi],3/4]]],
			GateDuration -> Abs[\[Theta]]/(freqcrot[n]) 
		|>,
		Subscript[CRy, e_,n_][\[Theta]_] /; (e==0 && n>0):> <|
			(* Its noisy form depolarises the control and target qubits *)
			NoisyForm -> Circuit[Subscript[CRy, e,n][\[Theta]]Subscript[Depol,e,n][Min[ercrot[n][[1]]Abs[\[Theta]]/\[Pi],15/16]]Subscript[Deph,e,n][Min[ercrot[n][[2]]Abs[\[Theta]]/\[Pi],3/4]]],
			GateDuration -> Abs[\[Theta]]/(freqcrot[n])    
		|>		
	}	
	,
	(* Declare that \[CapitalDelta]t will refer to the duration of the current gate/channel. *)
	DurationSymbol -> \[CapitalDelta]t, 

(****)
(* PASSIVE NOISE *)
(****)
	(* The passive noise on qubits when NOT being operated upon.  *)
(* Note 'globalField' is the unwanted [positive or negative] field offset for the present circuit; *)
(* this should be set on a per-run basis as in examples below. Could be augmented with e.g. a drifting function of t *)
(** other potential passive noise: Global magnetic field fluctuation due to ^13C bath 
Causes dephasing on the electron in 4 \[Mu]s
External field is effectively perfect
**)

		Qubits -> {
		q_ :> <|
				PassiveNoise -> pn[q,\[CapitalDelta]t]
			|>	
	}
|>

]
]

(**************************** RYDBERG_HUB_HARVARD *****************************)
(**legitimate shift move **) 
legShift[q_,v_,atomlocs_]:=Module[{qlocs=atomlocs},
	qlocs[#]+=v&/@Flatten[{q}];
	Length@qlocs === Length@DeleteDuplicates@Values@qlocs
];

asymBitFlip::usage="Asymmetric bit-flip error. Realised with amplitude damping and symmetric bitflip";
Subscript[asymBitFlip, q_][p01_,p10_]:=Module[{pbf=Min[p01,p10], pmax=Max[p01,p10],edamp,x},
	(*Here, we assume damping and bitflip are independent events. 
	P(damp or bf)=P(damp)+P(bf)-P(damp)P(bf)*)
	edamp=First[x/.Solve[x+pbf-x*pbf==pmax,{x}]];

	Which[
		(* more states are flipped to 1*)
		p01>p10,
		Sequence@@{Subscript[X, q],Subscript[Damp, q][edamp],Subscript[X, q],Subscript[Kraus, q][bitFlip1[1-pbf]]}
		,
		p01<p10,
		(*more states are flipped to 0*)
		Sequence@@{Subscript[Damp, q][edamp],Subscript[Kraus, q][bitFlip1[1-pbf]]}
		,
		True,
		(*equals*)
		Subscript[Kraus, q][bitFlip1[1-pbf]]
	]
]	
			
							
RydbergHub[OptionsPattern[]]:=With[
{
	qubitsnum=OptionValue@qubitsNum,
	atomlocations=OptionValue@AtomLocations,
	t2=OptionValue@T2,
	vaclifetime=OptionValue@VacLifeTime,
	rabifreq=OptionValue@RabiFreq,
	unitlattice=OptionValue@UnitLattice,
	blockaderad=OptionValue@BlockadeRadius,
	probleakinit=OptionValue@ProbLeakInit,
	probbfrot01=OptionValue[ProbBFRot][01],
	probbfrot10=OptionValue[ProbBFRot][10],
	durinit=OptionValue@DurInit,
	fidmeas=OptionValue@FidMeas,
	durmeas=OptionValue@DurMeas,
	problossmeas=OptionValue@ProbLossMeas,
	probleakcz=OptionValue@ProbLeakCZ,
	probleakmove=OptionValue@ProbLeakMove,
	qubits=Keys@OptionValue@AtomLocations
},

	(**** assertions ****)
	Catch@If[Length@atomlocations!=qubitsnum, Throw@Message[qubitsNum::error,"Missing or extra qubits in AtomLocations"]];
	Catch@If[\[Not](0<=probleakinit<=1), Throw@Message[ProbLeakInit::error,"Needs value within [0,1]"]];
	Catch@If[\[Not](0<=problossmeas<=1), Throw@Message[ProbLossMeas::error,"Needs value within [0,1]"]];
	Catch@If[\[Not](0<=probleakcz<=1), Throw@Message[ProbLeakCZ::error,"Needs value within [0,1]"]];
	Catch@If[\[Not](0<=probleakmove<=1), Throw@Message[ProbLeakMove::error,"Needs value within [0,1]"]];

Module[{\[CapitalDelta]t, lossatoms, lossatomsprob, globaltime, stdpn, t1, atomlocs, distloc,blockadecheck},
	atomlocs=atomlocations;
	(* the atoms that are lost to the environment *)
	lossatoms=<|Table[k->False,{k,Keys@atomlocations}]|>;
	lossatomsprob=<|Table[k->0,{k,Keys@atomlocations}]|>;
	t1=vaclifetime/qubitsnum;
	
	(** standard passive noise **) 
	stdpn[q_,dur_]:={Subscript[Depol, q][0.75(1-E^(-dur/t1))],Subscript[Deph, q][0.5(1-E^(-dur/t2))]};					
	(** coordinate distance measure *)
	distloc[q1_,q2_]:=Norm[atomlocs[q1]-atomlocs[q2],2]*unitlattice;
	
	(** legitimate multi-qubit gates**) 
	blockadecheck[q_List]:=And@@((distloc@@#<= blockaderad)&/@Subsets[Flatten[q],{2}]);

<|
	DeviceDescription -> ToString[qubitsnum]<>" Rydberg qubits in a 3D lattice.",
	NumAccessibleQubits -> qubitsnum,
	NumTotalQubits ->qubitsnum,

	(** custom keys to access internal variables **)
	LossAtomsProbability:>lossatomsprob,
	LossAtoms:>lossatoms,
	AtomLocations:>atomlocs,
	BlockadeRadius->blockaderad,
	UnitLattice->unitlattice,
	
	(*(* re-initialized when invoking InsertCircuitNoise *)
	InitVariables->Function[
		atomlocs=atomlocations;
		lossatoms=<|Table[k->False,{k,Keys@atomlocations}]|>;
		lossatomsprob=<|Table[k->0,{k,Keys@atomlocations}]|>;
	],*)

	(* Aliases are useful for declaring custom operators. At this stage the specification is noise-free (but see later) *)
	Aliases -> {
		Subscript[Init, q_Integer]:> {}
		,
		Subscript[SRot, q_Integer][\[Phi]_,\[CapitalDelta]_,tg_]:>Circuit[Subscript[UG, q][\[Phi],\[CapitalDelta],tg,rabifreq]]
		,
		Subscript[CZ, p_Integer,q_Integer][\[Phi]_]:>Circuit[Subscript[U, p,q][{{1,0,0,0},{0,E^(I \[Phi]),0,0},{0,0,E^(I \[Phi]),0},{0,0,0,E^(I(2 \[Phi]-\[Pi]))}}]]
		,
		Subscript[SWAPLoc, q1_Integer,q2_Integer]:>{}
		,
		Subscript[Wait, q__][t_] :>{}
		,
		Subscript[ShiftLoc, q__][v_]:>{}
		,
		(* multi-control-singleZ *)
		Subscript[C, c_Integer][Subscript[Z, t__Integer][\[Theta]_]]:>Table[Subscript[C, c][Subscript[Z, targ]],{targ,{t}}]
		,
		Subscript[C, c_Integer][Subscript[Z, t__Integer]]:>Table[Subscript[C, c][Subscript[Z, targ]],{targ,{t}}]
		,
		(* single-control multi-Z*)
		Subscript[C, c__Integer][Subscript[Z, t_Integer][\[Theta]_]]:>{Subscript[C, c][Subscript[Z, t]]}
	},
	
	(* Global time: not yet used here *)
	TimeSymbol-> globaltime,
	
	(* gates rules *)
	Gates -> {
	Subscript[Init, q_Integer] :> <|
	 (* Put the electron back to the atom and reset leak probability *)
		UpdateVariables-> Function[
						lossatoms[q]=False;
						lossatomsprob[q]=0],
		NoisyForm -> {Subscript[Damp, q][1],Subscript[KrausNonTP,q][{{{Sqrt[1-probleakinit],0},{0,1}}}]}, 
		GateDuration -> durinit
	|>
	,
	Subscript[Wait, q__][dur_]/;(Complement[Flatten@{q},qubits]==={}):> <|
	NoisyForm->Table[stdpn[j,dur],{j,Flatten@{q}}],
	GateDuration->dur
	|>
	,
	Subscript[M, q_Integer]:> <|
	UpdateVariables-> Function[
			lossatomsprob[q]=1-(1-lossatomsprob[q])*(1-problossmeas)],
	NoisyForm -> {Subscript[Depol, q][Min[1-fidmeas,3/4]],Subscript[M, q]},  
	GateDuration -> durmeas
	|>
	(** Single-qubit gates **)
	,
	Subscript[SRot, q_Integer][\[Phi]_,\[CapitalDelta]_,tg_]:><|
	NoisyForm->{Subscript[SRot, q][\[Phi],\[CapitalDelta],tg],Subscript[Deph, q][0.5(1-E^(-tg*0.5/t2))]},
	GateDuration-> tg
	|>
	,
	Subscript[H, q_Integer]:> <|
	NoisyForm->{Subscript[H, q],Subscript[Deph, q][0.5(1-E^(-0.5/(rabifreq*t2)))]},
	GateDuration-> \[Pi]/rabifreq
	|>
	,
	Subscript[Rx, q_Integer][\[Theta]_]:><|
	NoisyForm->{Subscript[Rx, q][\[Theta]],Subscript[asymBitFlip, q][probbfrot01,probbfrot10],Subscript[Deph, q][0.5(1-E^(-0.5*Abs[\[Theta]/\[Pi]]/(rabifreq*t2)))]},
	GateDuration-> Abs[\[Theta]]/rabifreq
	|>
	,
	Subscript[Ry, q_Integer][\[Theta]_]:><|
	NoisyForm->{Subscript[Ry, q][\[Theta]],Subscript[asymBitFlip, q][probbfrot01,probbfrot10],Subscript[Deph, q][0.5(1-E^(-0.5*Abs[\[Theta]/\[Pi]]/(rabifreq*t2)))]},
	GateDuration-> Abs[\[Theta]]/rabifreq
	|>
	,
	Subscript[Rz, q_Integer][\[Theta]_]:><|
	NoisyForm->{Subscript[Rz, q][\[Theta]],Subscript[asymBitFlip, q][probbfrot01,probbfrot10],Subscript[Deph, q][0.5(1-E^(-0.5*Abs[\[Theta]/\[Pi]]/(rabifreq*t2)))]},
	GateDuration-> Abs[\[Theta]]/rabifreq
	|>
	, 
	(** two-qubit gates **)
	Subscript[SWAPLoc, q1_Integer,q2_Integer]:> <|
	UpdateVariables-> Function[{atomlocs[q1],atomlocs[q2]}={atomlocs[q2],atomlocs[q1]};],
	NoisyForm-> Flatten@{stdpn[q1,4\[Pi]/rabifreq],stdpn[q2,4\[Pi]/rabifreq]},
	GateDuration->4\[Pi]/rabifreq
	|>
	,
	Subscript[ShiftLoc, q__][v_]/;legShift[Flatten@{q},v,atomlocs]:> <|
	UpdateVariables-> Function[atomlocs[#]+=v &/@Flatten[{q}];],
	NoisyForm-> Join[stdpn[#,4\[Pi]/rabifreq]&/@Flatten[{q}], Subscript[KrausNonTP,#][{{{1,0},{0,Sqrt[1-probleakmove]}}}]&/@Flatten[{q}]],
	GateDuration->4\[Pi]/rabifreq
	|>
	,
	Subscript[CZ, p_Integer,q_Integer][\[Phi]_]/;blockadecheck[{p,q}]:><|
	NoisyForm->{Subscript[CZ, p,q][\[Phi]],Subscript[KrausNonTP, p,q][{{{1,0,0,0},{0,Sqrt[1-probleakcz[01]],0,0},{0,0,Sqrt[1-probleakcz[01]],0},{0,0,0,Sqrt[1-probleakcz[11]]}}}]},
	GateDuration-> Abs[\[Phi]]/rabifreq
	|>
	,
	(* If argument is presented, it acts as gate duration per \[CapitalOmega]  *)
	Subscript[C, c_][Subscript[Z, t__]]/;blockadecheck[{c,t}]:><|
	NoisyForm->Join[Table[Subscript[C, c][Subscript[Z, targ]],{targ,{t}}], Subscript[KrausNonTP, #][{{{1,0},{0,Sqrt[1-probleakcz[11]]}}}]&/@Flatten@{c,t}],
	GateDuration->4\[Pi]/rabifreq
	|>
	,
	Subscript[C, c_][Subscript[Z, t__][\[Theta]_]]/;blockadecheck[{c,t}]:><|
	NoisyForm ->Join[Table[Subscript[C, c][Subscript[Z, targ]],{targ,{t}}], Subscript[KrausNonTP, #][{{{1,0},{0,Sqrt[1-probleakcz[11]]}}}]&/@Flatten@{c,t}],
	GateDuration-> Abs[\[Theta]]/rabifreq
	|>
	,
	Subscript[C, c__][Subscript[Z, t_]]/;blockadecheck[{c,t}]:> <|
	NoisyForm->Join[{Subscript[C, c][Subscript[Z, t]]},Subscript[KrausNonTP, #][{{{1,0},{0,Sqrt[1-probleakcz[11]]}}}]&/@Flatten@{c,t}],
	GateDuration->4\[Pi]/rabifreq
	|>
	,
	Subscript[C, c__][Subscript[Z, t_][\[Theta]_]]:> <|
	NoisyForm->Join[{Subscript[C, c][Subscript[Z, t]]},Subscript[KrausNonTP, #][{{{1,0},{0,Sqrt[1-probleakcz[11]]}}}]&/@Flatten@{c,t}],
	GateDuration->Abs[\[Theta]]/rabifreq
	|>
	}
	,
	DurationSymbol -> \[CapitalDelta]t,
	(* PASSIVE NOISE *)
		Qubits -> { q_ :> <| PassiveNoise -> stdpn[q,\[CapitalDelta]t]|>}
|>
	]
]
SetAttributes[PlotAtoms,HoldFirst]
PlotAtoms[rydbergdev_,opt:OptionsPattern[{PlotAtoms,Graphics}]]:=With[
{
	qulocs=rydbergdev[AtomLocations],
	blrad=rydbergdev[BlockadeRadius],
	unit=rydbergdev[UnitLattice],
	blockade=OptionValue[ShowBlockade],
	showloss=OptionValue[ShowLossAtoms],
	availqubits=KeyTake[rydbergdev[AtomLocations],Keys@Select[rydbergdev[LossAtoms],#===False&]],
	lossqubits=KeyTake[rydbergdev[AtomLocations],Keys@Select[rydbergdev[LossAtoms],#===True&]],
	style={ImageSize->Medium,Frame->True,Axes->True}(*,
	fontsize=If[(OptionValue[ImageSize]===Automatic)||NumberQ[OptionValue[ImageSize]],Medium,OptionValue[ImageSize]]*)
},
Which[
	2===Length@First@Values@qulocs,
	Show[
		Sequence@@Table[Graphics[{Cyan,Opacity[0.15],EdgeForm[Directive[Dashed,Orange]],Disk[unit*qulocs[b],blrad]}],{b,blockade}],
		Sequence@@Table[Graphics[{Red,Disk[unit*v,0.15]}],{v,Values@availqubits}],
		Sequence@@Table[Graphics[Text[k,({0.1,0.1}+qulocs[k])*unit]],{k,Keys@availqubits}],
		If[showloss,Sequence@@Table[Graphics[{Gray,Disk[unit*v,0.15]}],{v,Values@lossqubits}],Sequence@@{}],
		If[showloss,Sequence@@Table[Graphics[Text[k,unit*({0.15,0.15}+qulocs[k])]],{k,Keys@lossqubits}],Sequence@@{}],
			Evaluate@FilterRules[{opt}, Options[Graphics]],Sequence@@style,AxesLabel->{"x","y"}
	]
	,
	3===Length@First@Values@qulocs,
	Show[
		Sequence@@Table[Graphics3D[{Red,Sphere[v*unit,0.1]}],{v,Values@availqubits}],
		Sequence@@Table[Graphics3D[{Text[k,unit*({0.15,0.15,0.15}+qulocs[k])]}],{k,Keys@availqubits}],
		Sequence@@Table[Graphics3D[{Cyan,Opacity[0.1],Sphere[unit*qulocs[b],blrad]}],{b,blockade}],
		If[showloss,Sequence@@Table[Graphics3D[{Gray,Sphere[unit*v,0.1]}],{v,Values@lossqubits}],{}],
		If[showloss,Sequence@@Table[Graphics3D[Text[k,{0.15,0.15,0.15}+unit*qulocs[k]]],{k,Keys@lossqubits}],Sequence@@{}],
		Evaluate@FilterRules[{opt}, Options[Graphics]],Sequence@@style,AxesLabel->{"x","y","z"}
	]
	,
	True
	,
	Throw[Message[PlotAtoms::error,"Only accepts 2D- or 3D-coordinates."]]
	]
]
(****************** ENDOF RYDBERG_HUB*********************)

(***** SILICON_DELFT *****)
SiliconDelft[OptionsPattern[]]:=With[
{
	(*validate and format parameter specification*)
	qubitsnum=Catch@validate[OptionValue@qubitsNum,IntegerQ,qubitsNum,"not an integer"],
	(*Fractions*)
	efsinglexy=Catch@validate[OptionValue@EFSingleXY,(Total[#]==1||Total[#]==0)&,EFSingleXY,"not a fraction with total 1 "],
	efcz=Catch@validate[OptionValue@EFCZ,(Total[#]==1||Total[#]==0)&,EFCZ,"not a fraction with total 1 "],
	
	(*Number as average or association to specify each*)
	t1=Catch@validate[OptionValue@T1,checkAss[#,OptionValue@qubitsNum]&,T1,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
	t2=Catch@validate[OptionValue@T2,checkAss[#,OptionValue@qubitsNum]&,T2,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
	rabifreq=Catch@validate[OptionValue@RabiFreq,checkAss[#,OptionValue@qubitsNum]&,RabiFreq,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
	qubitfreq=Catch@validate[OptionValue@QubitFreq,checkAss[#,OptionValue@qubitsNum]&,QubitFreq,numass@OptionValue@qubitsNum,num2Ass[#,OptionValue@qubitsNum]&],
	freqcz=Catch@validate[OptionValue@FreqCZ,checkAss[#,-1+OptionValue@qubitsNum]&,FreqCZ,numass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
	
	(*Number as average or association to specify each fidelity*)
	fidcz=Catch@validate[OptionValue@FidCZ,checkAss[#,-1+OptionValue@qubitsNum,0<=#<=1&]&,FidCZ,fidass[-1+OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
	fidsinglexy=Catch@validate[OptionValue@FidSingleXY,checkAss[#,OptionValue@qubitsNum,0<=#<=1&]&,FidSingleXY,fidass[OptionValue@qubitsNum],num2Ass[#,-1+OptionValue@qubitsNum]&],
	
	(*single things*)
	fidread=Catch@validate[OptionValue@FidRead,0<=#<=1&,FidRead,"invalid fidelity"],
	durread=Catch@validate[OptionValue@DurRead,NumberQ,DurMeas,"invalid duration"],
	(*repeatread=Catch@validate[OptionValue@RepeatRead,IntegerQ,RepeatRead,"not an integer"],*)
	(* assoc or boolean *)
	exchangerotoff=Catch@validate[OptionValue@ExchangeRotOff,Or[AssociationQ@#,#===False]&,ExchangeRotOff,"Set to association or False"],
	
	(* True/False*)
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
{\[CapitalDelta]T, miseq,initf,measf, passivenoisecirc, offresrabi, stdpn, exczon,durinit,sroterr,g2=False,ndeph,ndepol},
	(* Normalise the numbers to be within the correct range of error parameters*)		
	ndeph[num_]:=Min[num,0.5];
	ndepol[num_]:=Min[num,0.75];	
	
	stdpn[q_,dur_]:=If[stdpassivenoise,{Subscript[Deph, q][.5(1-Exp[-N[dur/t2[q],$MachinePrecision]])],Subscript[Depol, q][0.75(1-Exp[-dur/t1[q]])]},{}];
	passivenoisecirc[q_Integer,g2_,dur_]:=Flatten@{If[\[Not]g2&&q<qubitsnum-1 && AssociationQ@exchangerotoff,Subscript[C, q][Subscript[Rz, q+1][(dur/\[Pi])*exchangerotoff[q]]],{}],stdpn[q,dur]};
	
	(*single rotation noise *)
	offresrabi[q_,\[Theta]_]:=If[offresonantrabi,Table[Subscript[U, j][OffResRabiOsc[rabifreq[q],qubitfreq[j]-qubitfreq[q],Abs[\[Theta]]]],{j,Delete[qubits,q+1]}],{}];
	
	(*Exchange rotation C-Rz[j] interaction when CZ gate on*)
	exczon[targ_]:=If[ListQ@exchangeroton,Subscript[C, #-1][Subscript[Rz, #][exchangeroton[[targ,#]]]]&/@Delete[Range[qubitsnum-1],targ],{}];
	
	(*Errors on single rotations*)
	sroterr[q_,\[Theta]_]:=Flatten@{Subscript[Depol, q][ndepol[er1xy[q][[1]]*Abs[\[Theta]/\[Pi]]]],Subscript[Deph, q][ndeph[er1xy[q][[2]]*Abs[\[Theta]/\[Pi]]]],offresrabi[q,\[Theta]]};
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
durinit[q__]:=durread*3+1.5*Total[Flatten@Table[rabifreq[i],{i, Complement[{q},{0,qubitsnum-1}]}]]; 
						
(* measurment has bitflip errors at the edge and single qubit errors in the middle *)		
measf[q__]:=Which[
	MemberQ[{q},0],(*start*)
	Flatten@{
	(*mimics 2 readout  *)
	{Subscript[Kraus, 0][bitFlip1[1-(1-fidread)^2]],Subscript[M, 0],Subscript[M, 1],Subscript[Damp, 0][1-(1-fidread)^2],Subscript[X, 0],Subscript[Damp, 1][1-(1-fidread)^2]},
	(* 2 rotation errors + 1 readout *)
	Table[{Subscript[Depol, i][er1xy[i][[1]]],Subscript[Deph, i][er1xy[i][[2]]],Subscript[M, i],Subscript[Damp, i][1-(1-fidread)]},{i,Complement[{q},{0,1}]}]
	}
	,
	MemberQ[{q},qubitsnum-1],(*end*)
	Flatten@{
	(*mimics 2 readout  *)
	{Subscript[Kraus, qubitsnum-1,qubitsnum-2][bitFlip2[1-(1-fidread)^2]],Subscript[M, qubitsnum-1],Subscript[M, qubitsnum-2],Subscript[Damp, qubitsnum-1][1-(1-fidread)^2],Subscript[X, qubitsnum-1],Subscript[Damp, qubitsnum-2][1-(1-fidread)^2]},
	(* 2 rotation errors + 1 readout *)
	Table[{Subscript[Depol, i][er1xy[i][[1]]],Subscript[Deph, i][er1xy[i][[2]]],Subscript[M, i],Subscript[Damp, i][1-(1-fidread)]},{i,Complement[{q},{qubitsnum-1,qubitsnum-2}]}]
	}
	,
	True,
	Throw[Message[M::error,"Error in Measurement operation"]]
];

	<|
	(*no hidden qubits/ancilla here *)
	DeviceType->"SiliconDelft",
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
			NoisyForm->Flatten@Table[passivenoisecirc[i,False,t],{i,Flatten@{q}}],
			GateDuration->t,
			UpdateVariables->Function[g2=False]
		|>,
	(* Measurements and initialisation *)
		Subscript[Meas, q__]/; miseq[q] :><|
			NoisyForm-> measf[q],
			GateDuration->durread,
			UpdateVariables->Function[g2=False]
		|>,
		Subscript[Init, q__]/; miseq[q] :><|
			NoisyForm-> initf[q],
			GateDuration-> durinit[q],
			UpdateVariables->Function[g2=False]
		|>,		
	(* Singles *)
		Subscript[Rx,q_][\[Theta]_]:><|
			NoisyForm->Flatten@{Subscript[Rx, q][\[Theta]],sroterr[q,\[Theta]]},
			GateDuration->Abs[\[Theta]]/(rabifreq[q]),
			UpdateVariables->Function[g2=False] 
		|>,
		Subscript[Ry,q_][\[Theta]_]:><|
			NoisyForm->Flatten@{Subscript[Ry, q][\[Theta]],sroterr[q,\[Theta]]},
			GateDuration->Abs[\[Theta]]/(rabifreq[q]),
			UpdateVariables->Function[g2=False]
		|>,
	(* Twos *)
		Subscript[C, p_][Subscript[Z, q_]]/; q-p===1  :><|
			(*The last bit undo the exchange in the passive noise *)
			NoisyForm->{Subscript[C, p][Subscript[Z, q]],Subscript[Depol, p,q][ercz[p][[1]]],Subscript[Deph, p,q][ercz[p][[2]]],Sequence@@exczon[q]}, 
			GateDuration->\[Pi]/freqcz[p],
			UpdateVariables->Function[g2=True]
		|>,
		Subscript[C, p_][Subscript[Ph, q_][\[Theta]_]]/; q-p===1  :><|
			(*The last bit undo the exchange in the passive noise *)
			NoisyForm->{Subscript[C, p][Subscript[Ph, q][\[Theta]]],Subscript[Depol, p,q][Min[ercz[p][[1]]Abs[\[Theta]/\[Pi]],15/16]],Subscript[Deph, p,q][Min[ercz[p][[2]]*Abs[\[Theta]/\[Pi]],3/4]],Sequence@@exczon[q]}, 
			GateDuration->Abs[\[Theta]]/freqcz[p],
			UpdateVariables->Function[g2=True]
		|>		
	},
	(* Declare that \[CapitalDelta]t will refer to the duration of the current gate/channel. *)
	DurationSymbol -> \[CapitalDelta]T,
	 
	(* Passive noise *)
	Qubits :> {
			q_/;(0<=q<qubitsnum) :> <|
			PassiveNoise ->passivenoisecirc[q,g2,\[CapitalDelta]T]
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
showIons[nodes_,connall_]:=Module[
{qulocs,i,colors},
	colors=rainbowcol[Length@nodes];
	qulocs=<|Table[node-><|Table[i=0;zone-><|Table[qubit->{i++,0},{qubit,nodes[node][zone]}]|>,{zone,Keys@nodes[node]}]|>,{node,Keys@nodes}]|>;
	i=0;
	Grid[Table[
	i++;{node,Sequence@@Table[
	drawZone[qulocs[node][zone],"zone "<>ToString@zone,colors[[i]],connall[node]],{zone,Keys@nodes[node]}]
	},{node,Keys@nodes}],Spacings->{0,1}]
]
rainbowcol[n_]:=Table[ColorData["CMYKColors",(i-1)/n],{i,n}]
drawZone[locszone_,label_,color_,conn_]:=Module[{v,edges},
edges=Select[Table[If[EdgeQ[conn,#\[UndirectedEdge]#2&[Sequence@@pair]],pair,False],{pair,Subsets[Keys@locszone,{2}]}],ListQ];
Show[
	Sequence@@Table[Graphics[{Dashed,Black,Line[{locszone[#],locszone[#2]}&[Sequence@@e]]}],{e,edges}],
	Graphics[Table[{color,Ball[loc,.25]},{loc,Values@locszone}]],
	Graphics[Table[Text[k,locszone[k],BaseStyle->"Output"],{k,Keys@locszone}]],
	ImageSize->{50*Max[1,Length@locszone],50}, Frame->True,FrameTicks->False,FrameStyle->Black,FrameLabel->label,FrameMargins->Tiny
	]
]
(* zone-related functions *)
getZone[q_,node_]:=(Association@Flatten@Table[Table[v->k,{v,node[k]}],{k,Keys@node}])[q]

(* 
Legitimate split moves:
	1) both ions initially stay in the same zone: zone 1-3.
	2) destination split must be different with the initial zone.
	3) the ions assigned to split must stay next to each other.
	4) there is no ions in between zones: the split is linear, no jumps.
*)
legSplit[nodes_,nodename_,q1_,q2_,zone_]:=With[{z1=getZone[q1,nodes[nodename]],z2=getZone[q2,nodes[nodename]],node=nodes[nodename]},
	If[z1!=z2,Return@False]; If[(z1===4)||(z2===4),Return@False];
	If[(z1===zone)||(z2===zone),Return@False];
	If[Abs[Position[node[z1],q1][[1,1]]- Position[node[z1],q2]][[1,1]]!=1,Return@False];
	If[Length@Flatten[Table[node[z],{z,Range[Min[z1,zone]+1,Max[z1,zone]-1]}]]>0,Return@False];
	True
]

(* Legitimate shuttling moves:
	1) The zone destination is different with the start zone
    2) Shuttle only works for the entire ions sitting in a zone
    3) There is no ions in between moves, shuttle is linear
*)
legShutl[nodes_,nodename_,zone_,qubits__]:=With[
{zstart=getZone[#,nodes[nodename]]&/@{qubits},node=nodes[nodename]},
If[\[Not]Equal@@zstart,Return@False];
If[nodes[nodename][zstart[[1]]]!={qubits},Return@False];
Length@Flatten[Table[node[z],{z,Range[Min[zone,zstart[[1]]]+1,Max[zone,zstart[[1]]]-1]}]]===0
]

(* Legitimate combine moves:
1) When the final zone unspecified, the ions must be sitting in the same zone
2) One can move to a higher zone or lower zone, zone destination must be 1-3. 
3) No ions in between, combine shuttle is linear
*)
legComb[nodes_,nodename_,q1_,q2_,zonedest_:None,connectivity_:None]:=Module[{zone,node=nodes[nodename],z1,z2,cond1,cond2,zstart,ps,pz,qs,qz,cond3},
	(* combine to the one of the zone *)
	z1=getZone[q1,node];
	z2=getZone[q2,node];
	
	(* Unspecified zone destination must be done within the same zone*)
	zone=If[zonedest===None, z1, zonedest];
	cond1=Or[z1===zone, z2===zone];
	
	(* [check again] *)
	If[\[Not]cond1, Return[False]];
	If[z1===zone,
		zstart=z2;qs=q2;qz=q1,
		zstart=z1;qs=q1;qz=q2
	];
	If[zone===4, Return[False]];	
	(* merge from start zone to zone destination *)
	ps=Position[node[zstart],qs][[1,1]];
	pz=Position[node[zone],qz][[1,1]];
	cond2=Which[
	(* move down to a higher zone *)
	zstart<zone,
		Length@node[zstart][[ps+1;;]]+Length@node[zone][[;;pz-1]](*+Length@EdgeList@NeighborhoodGraph[connectivity[nodename],qs]*),
	zstart>zone,
	(*move up to a lower zone *)
		Length@node[zone][[pz+1;;]]+Length@node[zstart][[;;ps-1]](*+Length@EdgeList@NeighborhoodGraph[connectivity[nodename],qs]*),
	zstart===zone,
	(* doesn't move *)
	Length@node[zone][[Min[ps,pz]+1;;Max[ps,pz]-1]]
	];
	
	If[cond2!=0,Return[False]];

	(* no qubits in between zones *)
	cond3=If[zone===zstart,
	True,
	Length@Flatten[Table[node[z],{z,Range[Min[zstart,zone]+1,Max[zstart,zone]-1]}]]===0];
	If[\[Not]cond3, Return[False]];
	True
]

(** Split move **)
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

(* move qubit q to a zone*)
Subscript[shutl, q__][inodes_,node_,zone_]:=Module[{szone, nodes},
	nodes=inodes;
	szone=getZone[{q}[[1]],nodes[node]];
	If[zone>szone,
	nodes[node][zone]=Join[nodes[node][szone],nodes[node][zone]]
	,
	nodes[node][zone]=Join[nodes[node][zone],nodes[node][szone]]
	];
	nodes[node][szone]={};
	nodes
]

Subscript[comb, i_,j_][inodes_,nodename_,zone_]:=Module[{nodes,node,z1,z2,ps,pz,sq,zstart,qs,qz},
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
	bfprob=OptionValue@BFProb,
	durinit=OptionValue@DurInit,
	durmove=OptionValue@DurMove,
	durread=OptionValue@DurRead,
	efent=OptionValue@EFEnt,
	fidinit=OptionValue@FidInit,
	freqcz=OptionValue@FreqCZ,
	freqent=OptionValue@FreqEnt,
	fident=OptionValue@FidEnt,
	rabifreq=OptionValue@RabiFreq,
	scatprob=OptionValue@ScatProb,
	t1=OptionValue@T1,
	t2s=OptionValue@T2s,
	stdpn=OptionValue@StdPassiveNoise,
	
	(*probability error of depolarising and dephasing noise *)
	er1xy=Association[#->fid2DepolDeph[OptionValue[FidSingleXY][#],OptionValue[EFSingleXY][#],1,FidSingleXY,True]&/@Keys[OptionValue@Nodes]], 
	ercz=Association[#->fid2DepolDeph[OptionValue[FidCZ][#],OptionValue[EFCZ][#],2,FidCZ,True]&/@Keys[OptionValue@Nodes]],
	
	(* convert the entanglement fidelity to the gate fidelity first here *)
	erent=entfid2DepolDeph[OptionValue[FidEnt],OptionValue[EFEnt],FidEnt]
	},
	
	Module[
	{\[CapitalDelta]t, nodes, qmap, qnum,checkpsd,checkread,checklog,checkrent,gnoise,entnoise,swaploc,connectivity,passivenoise,scatnoise,conninit,nodesinit,stdp},
	
	{nodesinit,qmap,qnum}=createNodes[initnodes];
	(* empty list of connectivitiy of two-qubit gates generated by Comb gate. stored as graphs *)
	conninit := Association[Table[n->With[{vex=nodes[n][1]},Graph[Table[vex[[i]]\[UndirectedEdge]vex[[i+1]],{i,-1+Length@vex}],VertexLabels->"Name"]],{n,Keys@nodes}]];

	nodes=nodesinit;
	connectivity=conninit;
	
	(**scattering on the neighborhood **) 
	scatnoise[qubit_, node_]:=With[{g=Graph@ReplaceList[getZone[qubit,nodes[node]],{p___,a_,b_,q___}:>a\[UndirectedEdge]b]}, 
		Sequence@@Table[Subscript[Depol, n][scatprob[node]],{n,AdjacencyList[g,qubit]}]]; 
	(** check zones **)
	(* prepare, store, detect *)
	checkpsd[q_,node_]:=MemberQ[{1,2,3},getZone[q,nodes[node]]];
	checkread[q_,node_]:=With[{zone=getZone[q,nodes[node]]},And[MemberQ[{1,2,3},zone], Length@nodes[node][zone]==1]];
	
	checkpsd[p_,q_,node_]:=And[MemberQ[{1,2,3},getZone[q,nodes[node]]],MemberQ[{1,2,3},getZone[p,nodes[node]]],EdgeQ[connectivity[node],p\[UndirectedEdge]q]];
	(*logic*)
	checklog[q_,node_]:=MemberQ[{2,3},getZone[q,nodes[node]]];
	checklog[p_,q_,node_]:=And[MemberQ[{2,3},getZone[q,nodes[node]]],getZone[q,nodes@node]===getZone[p,nodes@node],EdgeQ[connectivity[node],p\[UndirectedEdge]q]];
	
	(*remote entangle*)
	checkrent[q1_,node1_,q2_,node2_]:=With[{z1=getZone[q1,nodes[node1]],z2=getZone[q2,nodes[node2]]},
											And[4===z1,4===z2,Length[nodes[node1][z1]]==1,Length[nodes[node2][z2]]==1]
									   ];
	
	(**gate noise **)
	gnoise[q_,node_,\[Theta]_]:=Sequence@@{Subscript[Depol, qmap[node][q]][Min[er1xy[node][[1]]*Abs[\[Theta]/\[Pi]],0.75]],Subscript[Deph, qmap[node][q]][Min[er1xy[node][[2]]*Abs[\[Theta]/\[Pi]],0.5]]};
	gnoise[p_,q_,node_,\[Theta]_]:=Sequence@@{Subscript[Depol, qmap[node][p],qmap[node][q]][Min[ercz[node][[1]]*Abs[\[Theta]/\[Pi]],15/16]],Subscript[Deph, qmap[node][p],qmap[node][q]][Min[ercz[node][[2]]*Abs[\[Theta]/\[Pi]],3/4]]};
	(* entanglement noise *)
	entnoise[q1_,q2_,node1_,node2_]:=Sequence@@{Subscript[Depol, qmap[node1][q1],qmap[node2][q2]][erent[[1]]],Subscript[Deph, qmap[node1][q1],qmap[node2][q2]][erent[[2]]]}; 
	
	(** extra operations **)
	swaploc[p_,q_,node_]:=Module[{z=getZone[q,nodes[node]],pos,lst},
		lst=nodes[node][z];
		pos=Flatten@{Position[lst,p],Position[lst,q]};
		lst[[pos]]=lst[[Reverse@pos]];
		nodes[node][z]=lst;
		nodes
	];
	
	(** get a list of idling qubits **)
	(* passive noise: T1 exp decay and T2* with gaussian decay *)
	stdp[node_,t_,q_]:={Subscript[Damp, qmap[node][q]][.75(1-E^(-t/t1[node]))],Subscript[Deph, qmap[node][q]][.5(1-Exp[-(t/t2s[node])^2])]};
	passivenoise[node_,t_,q__]:=Sequence@@If[stdpn,Flatten@Table[stdp[node,t,i],{i,Complement[Flatten@Values[nodes[node]],{q}]}],{}];
	passivenoise[node_,t_]:=Sequence@@If[stdpn,Flatten@{stdp[node,t,#]&/@Flatten[Values[nodes[node]]]},{}];
<|
	(*no hidden qubits/ancilla here *)
	DeviceType->"TrappedIonOxford",
	DeviceDescription -> StringForm["Trapped ion device Oxford style with ``. nodes",Length@nodes],
	NumAccessibleQubits -> qnum,
	NumTotalQubits -> qnum,
	Nodes:>nodes,
	QMap:>qmap,
	ShowNodes:>showIons[nodes,connectivity],
	
	(* re-initialized when invoking InsertCircuitNoise *)
	(*InitVariables->Function[
		nodes=nodesinit;
		connectivity=conninit;	
	],*)
	
	(* Init, Read, Rx, Ry, C[Z], Ent[node1,node2], SWAP, Splz, Comb  *)
	Aliases -> {
		Subscript[Wait, q__][node_,t_] :> Sequence@@{},
		Subscript[Init, q_][fid_]:> Subscript[Damp, q][fid],
		Subscript[Splz, i_,j_]:> Sequence@@{},
		Subscript[Comb, i_,j_]:>Sequence@@{},
		Subscript[Read, q_]:>Subscript[M, q],
		Subscript[Ent, p_,q_]:>Sequence@@{Subscript[Damp, p][1],Subscript[Damp, q][1],Subscript[X, q],Subscript[H, p],Subscript[C, p][Subscript[X, q]]},
		Subscript[SWAPLoc, i_,j_][t_]:>Sequence@@{},
		Subscript[Shutl, q__]:>Sequence@@{}
		}
	,	
Gates ->{
	Subscript[Rx, q_][node_,\[Theta]_]/;checklog[q,node]:> With[{dur=Abs[\[Theta]]/rabifreq[node]}, 
	<|
		NoisyForm->{Subscript[Rx, qmap[node][q]][\[Theta]],gnoise[q,node,\[Theta]],passivenoise[node,dur,q]},  
		GateDuration->dur
	|> ],
	Subscript[Ry, q_][node_,\[Theta]_]/;checklog[q,node]:>With[{dur=Abs[\[Theta]]/rabifreq[node]},
	<|
		NoisyForm->{Subscript[Ry, qmap[node][q]][\[Theta]],gnoise[q,node,\[Theta]],passivenoise[node,dur,q]},
		GateDuration->dur
	|> ],
	(* effectively virtual, noiseless *)
	Subscript[Rz, q_][node_,\[Theta]_] :>
	<|
		NoisyForm->{Subscript[Rz, qmap[node][q]][\[Theta]]},
		GateDuration->0.
	|>
	,
	Subscript[CZ, i_,j_][node_]/;checklog[i,j,node]:> With[{dur=\[Pi]/freqcz[node]},
	<|
		NoisyForm->{Subscript[C, qmap[node][i]][Subscript[Z, qmap[node][j]]],gnoise[i,j,node,\[Pi]],passivenoise[node,dur,i,j]},
		GateDuration->dur
	|>	],
	Subscript[Ent, q1_,q2_][node1_,node2_]/;checkrent[q1,node1,q2,node2]:>With[{dur=1/freqent},
	<|
		NoisyForm-> {Subscript[Ent, qmap[node1][q1],qmap[node2][q2]],entnoise[q1,q2,node1,node2],passivenoise[node1,dur,q1],passivenoise[node2,dur,q2]},
		GateDuration->dur
	|>   ],
	Subscript[Read, q_][node_]/; checkread[q,node] :><|
		NoisyForm->{Subscript[Kraus, qmap[node][q]][bitFlip1[bfprob@node]],Subscript[Read, qmap[node][q]],passivenoise[node,durread[node],q]}, 
		GateDuration->durread[node]
	|>,
	Subscript[Init, q_][node_]/; checkpsd[q,node] :><|
		NoisyForm->{Subscript[Init, qmap[node][q]][fidinit[node]],scatnoise[q,node],passivenoise[node,durinit[node],q]}, 
		GateDuration->durinit[node]
	|>,
	Subscript[Wait, q__][node_,t_]:><|
		NoisyForm->{passivenoise[node,t]},
		GateDuration->t
	|>,
	(** physically moving operations **)
	Subscript[SWAPLoc, i_,j_][node_]/;checkpsd[i,j,node] :><|
		NoisyForm->{Subscript[SWAPLoc, qmap[node][i],qmap[node][j]][durmove[node]],passivenoise[node,durmove[node][SWAPLoc]]}, 
		GateDuration->durmove[node][SWAPLoc],
		UpdateVariables->Function[nodes=swaploc[i,j,node]]
	|>,
	
	Subscript[Shutl, q__][node_, zone_]/;legShutl[nodes,node,zone,q]:><|
		NoisyForm->{Subscript[Shutl, Sequence[qmap[node]/@{q}]],passivenoise[node,durmove[node][Splz]]},
		GateDuration->durmove[node][Shutl],
		UpdateVariables->Function[
					nodes=Subscript[shutl, q][nodes,node,zone]
						]
	|>,
	
	Subscript[Splz, i_,j_][node_,zone_]/;legSplit[nodes,node,i,j,zone] :><|
		NoisyForm-> Flatten@{Subscript[Splz, qmap[node][i],qmap[node][j]],passivenoise[node,durmove[node][Splz]]},
		GateDuration->durmove[node][Splz],
		UpdateVariables->Function[
		nodes=Subscript[splitZ, i,j][nodes,node,zone];
		If[EdgeQ[connectivity[node],i\[UndirectedEdge]j],connectivity[node]=EdgeDelete[connectivity[node],i\[UndirectedEdge]j]];
		]
	|>,
	Subscript[Comb, i_,j_][node_,zone_]/;legComb[nodes,node,i,j,zone,connectivity]:><|
		NoisyForm->Flatten@{Subscript[Comb, qmap[node][i],qmap[node][j]],passivenoise[node,durmove[node][Comb]]},
		GateDuration->durmove[node][Comb],
		UpdateVariables->Function[
		nodes=Subscript[comb, i,j][nodes,node,zone];
		If[\[Not]EdgeQ[connectivity[node],i\[UndirectedEdge]j],
		connectivity[node]=EdgeAdd[connectivity[node],i\[UndirectedEdge]j]]
		]	
	|>,
	(* combine within the same zone *)
		Subscript[Comb, i_,j_][node_]/;legComb[nodes,node,i,j]:><|
		NoisyForm->Flatten@{Subscript[Comb, qmap[node][i],qmap[node][j]],passivenoise[node,durmove[node][Comb]]},
		GateDuration->durmove[node][Comb],
		UpdateVariables->Function[
		nodes=Subscript[comb, i,j][nodes,node,getZone[i,nodes[node]]];
		If[\[Not]EdgeQ[connectivity[node],i\[UndirectedEdge]j],
		connectivity[node]=EdgeAdd[connectivity[node],i\[UndirectedEdge]j]]
		] 
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
(*******  ENDOF_TRAPPED_IONS_OXFORD   *****)

(******* CIRCUIT_CONSTRUCTION ***********)
SetAttributes[CircRydbergHub, HoldAll]
blockadeParallel::usage="BlockadeParallel[g1,g2,device]. Check if gates g1 and g2 can be done concurrently. Some gates, such as Wait and Init has no parallel restriction.";
SetAttributes[blockadeParallel,HoldAll]
blockadeParallel[gate1_,gate2_,nadevice_]:=Module[{g1,g2,idx1,idx2,parallelgates,serialgates,distloc,atomlocs,unitlattice},
	atomlocs=nadevice[AtomLocations];
	unitlattice=nadevice[UnitLattice];
	distloc[q1_,q2_]:=Norm[atomlocs[q1]-atomlocs[q2],2]*unitlattice;
	serialgates={SWAPLoc,ShiftLoc,Wait};
	parallelgates={Init,SRot,H,Rx,Ry,Rz};
	{g1,idx1}=gateIndex[gate1];
	{g2,idx2}=gateIndex[gate2];
	Which[MemberQ[parallelgates,(g1|g2)],
	True,
	MemberQ[serialgates,(g1|g2)],
	False,
	True,
	And@@(distloc[Sequence@@#,nadevice]> nadevice[BlockadeRadius]&/@Tuples[{idx1,idx2}])]
]
CircRydbergHub[circuit_, device_, OptionsPattern[]]:=Module[
	{
	parallel=OptionValue[Parallel],
	circ=circuit,
	newcirc={},
	circols,idxcol,incol,idx1,idx2,g1,g2
	},
	
	If[
	(** parallel mode **)
	parallel,
		While[Length@circ>0 ,
			circols=GetCircuitColumns[circ];
			(* update equivalent ordering of the circuit *)
			circ=Flatten@circols;
			(* get indices partitioned wrt circols *)
			idxcol=TakeList[Range[Length@circ],Length@#&/@circols][[1]];
			
			(* eliminate non-legitimate gates of the first column *)
			AppendTo[newcirc,{}];
			incol=<| #->True & /@idxcol |>;
			
			Table[
			If[incol[i1],
				(* add gate i1 and eliminate the rest *)
				AppendTo[newcirc[[-1]],circ[[i1]]];
				Table[
				If[incol[[i2]],
				incol[[i2]]=blockadeParallel[circ[[i1]],circ[[i2]],device];
				]
				,{i2,Complement[idxcol,{i1}]}];
			]
			,{i1,idxcol}];
			(* update circuit *)
			circ=Delete[circ,{#}&/@Keys@Select[incol,#&]];
		]
	,
	(** serial mode **)
	While[Length@circ>0 ,
		circols=GetCircuitColumns[circ];
		(* update equivalent ordering of the circuit *)
		circ=Flatten@circols;
		(* get indices partitioned wrt circols *)
		idxcol=TakeList[Range[Length@circ],Length@#&/@circols][[1]];
		
		(* eliminate non-legitimate gates of the first column *)
		AppendTo[newcirc,{}];
		incol=<| #->True & /@idxcol |>;
		
		Table[
		If[incol[i1],
		(* add gate i1 and eliminate the rest *)
		AppendTo[newcirc[[-1]],circ[[i1]]];
		Table[
		If[incol[[i2]],
			{g1,idx1}=gateIndex[circ[[i1]]];
			{g2,idx2}=gateIndex[circ[[i2]]];
		incol[[i2]]=MemberQ[{Init,Wait},(g1|g2)]
		];
		,{i2,Complement[idxcol,{i1}]}];
		]
		,{i1,idxcol}];
		(* update circuit *)
		circ=Delete[circ,{#}&/@Keys@Select[incol,#&]];
	];	
	];
	newcirc
]
gateToParties[gate_]:=gate/.{ Subscript[Splz, __][n_,_]:>{n},
						Subscript[Comb, __][n_,___]:>{n},
						Subscript[Shutl, __][n_,_]:>{n},
						Subscript[Init, _][n_]:>{n},
						Subscript[Read, _][n_]:>{n},
						Subscript[Ent, __][n1_,n2_]:>{n1,n2},
						Subscript[SWAPLoc, __][n_]:>{n},
						Subscript[Rx, _][n_,_]:>{n},
						Subscript[Ry, _][n_,_]:>{n},
						Subscript[CZ, __][n_]:>{n}}
SetAttributes[CircTrappedIons,HoldAll]
CircTrappedIons[circuit_,device_,OptionsPattern[]]:=Module[
	{ccirc,i,entp,parties,placed,sidx=1,
	sorted,modify,
	parallel=OptionValue[Parallel],
	mapq=OptionValue[MapQubits],
	qmap=device[QMap],
	template=<|Table[n->None,{n,Keys@device[Nodes]}]|>
	},
	sorted={};
	modify[ass_Association,node_,val_]:=Module[{a=ass},a[node]=val;a];
	
	AppendTo[sorted,template];
	Table[
		parties=gateToParties[gate];
		placed=False;
		For[i=sidx, i<=Length@sorted, i++,
			If[
			(*find the right place *)
			And@@(sorted[[i]][#]==None&/@parties),
				placed=i;
				Break[]
			];
		];
		(* the right place is not found *)
		If[placed==False,
			AppendTo[sorted,template];
			placed=Length@sorted
		];			
		sorted[[placed]]=modify[sorted[[placed]],First@parties,gate];
		If[Length@parties==2, 
			sorted[[placed]]=modify[sorted[[placed]],Last@parties,True]
		];
		(*update the start search *)
		If[\[Not]MemberQ[Values[sorted[[placed]]],None],
			sidx=placed;
		]
		
	,{gate,Flatten@circuit}];
	ccirc=Table[DeleteCases[DeleteCases[Values@c,None],True],{c,sorted}];
	If[mapq,
	ccirc/.{
	Subscript[Shutl, qs__][n_,z_]:>Subscript[Shutl, Sequence@@Table[qmap[n][q],{q,{qs}}]][n,z],
	Subscript[Splz, q1_,q2_][n_,z_]:>Subscript[Splz, qmap[n][q1],qmap[n][q2]],
	Subscript[Comb, q1_,q2_][n_,z_]:>Subscript[Comb, qmap[n][q1],qmap[n][q2]],
	Subscript[Comb, q1_,q2_][n_]:>Subscript[Comb, qmap[n][q1],qmap[n][q2]],
	Subscript[Ent, q1_,q2_][n1_,n2_]:>Subscript[Ent, qmap[n1][q1],qmap[n2][q2]],
	Subscript[CZ, q1_,q2_][n_]:>Subscript[C, qmap[n][q1]][Subscript[Z, qmap[n][q2]]],
	Subscript[Rx, q_][n_,\[Theta]_]:>Subscript[Rx, qmap[n][q]][\[Theta]],
	Subscript[Ry, q_][n_,\[Theta]_]:>Subscript[Ry, qmap[n][q]][\[Theta]],
	Subscript[Rz, q_][n_,\[Theta]_]:>Subscript[Rz, qmap[n][q]][\[Theta]],
	Subscript[SWAPLoc, q1_,q2_][n_]:>Subscript[SWAPLoc, qmap[n][q1],qmap[n][q2]],
	Subscript[Init, q_][n_]:>Subscript[Init, qmap[n][q]],
	Subscript[Read, q_][n_]:>Subscript[M, qmap[n][q]]
	},
	ccirc
	]	
]

CircSiliconDelft[circ_, device_,OptionsPattern[]]:=Module[
{parallel=OptionValue[Parallel]},
Which[
	parallel===False,
	List/@Flatten[circ],
	True,
	Throw[Message[CircSiliconDelft::error,"unknown/unset Parallel value:"<>ToString[parallel]]]
	]
]

Serialize[circ_List]:=List/@Flatten[circ]
	
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
GenerateOptionTable[options_,columnwidth_:{"3cm","3cm","10cm"}]:=Module[{vars, header,content,contenf,contenl,fvalues,finfo},
vars=Keys@options;
header=StringForm["\\begin{tabular}{p{``}p{``}p{``}}\n\\toprule",Sequence@@columnwidth];
fvalues=List@@#&/@Values@options;
finfo=ToString@Information[#]&/@Keys[options];
content=Transpose@{Keys@options,fvalues,finfo};
contenl=Join[{"\\textbf{Variable}&\\textbf{Value}&\\textbf{Description}"},StringRiffle[#,"&"]&/@content];
StringRiffle[Join[{header},contenl,{"\\bottomrule"}],"\\\\ \n"]<>"\n\\end{tabular}"
]
SetAttributes[symbolNameJoin, HoldAll];
symbolNameJoin[symbols__Symbol] := Symbol @ Apply[
  StringJoin,
  Map[
   Function[s, ToString[Unevaluated[s]], HoldFirst],
   Hold[symbols]
  ]
];

gateIndex::usage="gateIndex[gate], returns {the gate, the indices}.";
gateIndex[gate_]:=With[{
	gpattern={
	R[t_,Subscript[g1_, p_] Subscript[g2_, q_]]:> {symbolNameJoin[g1,g2], {p,q}}, 
	Subscript[C, p__][Subscript[g_, q__]]:>{symbolNameJoin[C,g], {Sequence@@Flatten@{p,q}}},
	Subscript[g_, q__]:> {g, {Sequence@@Flatten@{q}}},
	Subscript[g_, q__][_]:> {g, {Sequence@@Flatten@{q}}}
	}},
	gate/. gpattern
];

(* list of gates that are always parallel *)
parallelGates=<|
	"SiliconDelft"->{Wait},
	"SiliconHub"->{Wait},
	"SuperconductingFZJ"->{Wait},
	"SuperconductingHub"->{Wait},
	"TrappedIonOxford"->{Wait},
	"TrappedIonInnsbruck"->{Wait},
	"RydbergHub"->{Wait,SWAPLoc,ShiftLoc,Wait},
	"RydbergWisconsin"->{Wait},
	"NVCenterDelft"->{Wait},
	"NVCenterHub"->{Wait}
|>;

(*
Partial trace on n-qubit
1) reshape to the tensor:ConstantArray[2,2*n]
2) contract
3) reshape to matrix with dim2^mx2^mwhere m=(n-#contract)
*)
PartialTrace[\[Rho]_List, qubits___]:=ptrace[\[Rho],qubits]
PartialTrace[\[Rho]_Integer, qubits___]:=ptrace[GetQuregMatrix[\[Rho]],qubits]
ptrace[\[Rho]mat_List,qubits___]:=Module[{tmat,nq,pmat,nfin,pairs},
	If[Length@{qubits}<1,Return[\[Rho]mat]];
	nq=Log2@Length@\[Rho]mat;
	(*tensorize*)
	tmat=ArrayReshape[\[Rho]mat,ConstantArray[2,2*nq]];
	(* contraction pairs, beware the least significant bit convention!! *)
    pairs=Reverse@Table[{i,i+nq},{i,Range[nq]}];
	pmat=TensorContract[tmat,pairs[[Sequence@@#]]&/@(1+{qubits})];
	nfin=nq-Length@{qubits};
	ArrayReshape[pmat,{2^nfin,2^nfin}]
]

(*** fidelity of two density matrices (Tr[Sqrt[\[Rho]\[Sigma]]])^2 ***)
CalcFidelityDensityMatrices[\[Rho]_,\[Sigma]_]:=Re[Tr[MatrixPower[If[IntegerQ@\[Rho],GetQuregMatrix[\[Rho]],\[Rho]] . If[IntegerQ@\[Sigma],GetQuregMatrix[\[Sigma]],\[Sigma]],1/2]]^2]

End[];
EndPackage[];
Needs["VQD`ParameterDevices`"]

