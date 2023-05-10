(* ::Package:: *)

BeginPackage["VQD`"];

(*
NOTE,TODO: 
Pausing the process does not always give the correct behaviour for the following reasons:
	1) it doesn't remember previous circuit
	2) the time starts from 0. Solution: use globaltime
*)

(** Devices **)
(*silicon devices*)

(*
SiliconDelft::usage = "Returns device specification of a Silicon device based on the device built by the University of Delft.";
SiliconDelft2::usage = "Returns device specification of a Silicon device with twice error severity of SiliconDelft.";
*)
SiliconHub::usage = "Returns devices specification of a Silicon device based on the device built by the QCSHub.";

(*superconducting qubit devices*)
(*
SuperconductingFZJ::usage = "Returns device specification of a Superconducting qubit device based on the device built by Forschungzentrum Juelich.";
*)
SuperconductingHub::usage = "Returns device specification of a Superconducting qubit device based on the device built by the QCSHub.";

(*trapped ion devices*)
TrappedIonOxford::usage = "Returns device specification of a multi-nodes Trapped ions based on the device built by the Oxford/Hub.";
(*
TrappedIonInnsbruck::usage = "Returns device specification of a string of Trapped ions base on the device built by the University of Innsbruck.";
*)

(*rydberg quantum devices/neutral atoms*)
RydbergHub::usage = "Returns device specification of a Rydberg/Neutral Atom device based on the device built by the QCSHub.";
(*
RydbergWisconsin::usage = "Returns device specification of a Rydberg/Neutral Atom device based on the device built by the University of Wisconsin.";
*)

(* nuclear-vacancy center devices *)
NVCenterDelft::usage = "Returns device specification of a Nitrogen-Vacancy diamond center device based on the device built by the University of Delft.";
(*
NVCenterHub::usage = "Returns device specification of a Nitrogen-Vacancy diamond center device based on the device built by the QCSHub.";
*)
(* toy device *)
ToyDevice::usage = "Return a specification with simple standard model.";

(** General functions  **)
CalcFidelityDensityMatrices::usage = "CalcFidelityDensityMatrices[\[Rho],\[Sigma]] fidelity of two density matrices, \[Rho] and \[Sigma] can be density matrix of Quregs. Fidelity of two density matrices.";
PartialTrace::usage = "PartialTrace[qureg/density matrix, tracedoutqubits_List]. Return the partial trace as a matrix.";
RandomMixState::usage = "RandomMixState[nqubits, nsamples:None]. Return a random mixed quantum density state matrix.";
GenerateOptionTable::usage = "GenerateOptionTable[options,columnwidth:{3cm,3cm,10cm}]. Create summary of options in Latex.Very buggy.";

(** Information keys **)
OptionsUsed::usage = "Show all options used in a virtual device specification instance.";
DeviceType::usage = "The type of device. Normally, the name of the function that generates it.";

(***  Visualisations   ***)
(* trapped ions *)
DrawIons::usage = "Draw the current cunfiguration string of ions in ion traps.";
DrawIons::error = "`1`";
(* neutral atoms *)
PlotAtoms::usage = "PlotAtoms[rydberg_device]. Plot the neutral atoms. Set ShowBlockade->{atoms} to visualise the blockade radii. Set ShowLossAtoms->True, to show the atoms that are loss as well (grey color). It also receives options of Graphic function to style the visualisation.";
PlotAtoms::error = "`1`";
Options[PlotAtoms] = {ShowBlockade -> {}, ShowLossAtoms -> False};

ShowBlockade::usage = "List the qubits to draw the blockade radius.";
ShowLossAtoms::usage = "Set true to show the atoms lost into the evironment. This shows the last coordinate before being lost.";

(** Circuit arrangements **)
(* default arrangement according to devices *)
Options[CircTrappedIons] = {MapQubits -> True, Parallel -> False};
Options[CircSiliconDelft] = {Parallel -> False};
Options[CircRydbergHub] = {Parallel -> False};
CircTrappedIons::usage = "CircTrappedIons[circuit, device, MapQubits->True, Parallel->False]. Circuit arrangement according to the device. Note that Parallle->True is not available yet.";
CircSiliconDelft::usage = "CircSiliconDelft[circuit, device, Parallel->False]. Circuit arrangement according to the device. Note that Parallle->True is not available yet.";
CircRydbergHub::usage = "CircRydbergHub[circuit, device, Parallel->(False, True)]";
Serialize::usage = "Serialize circuit. Every quantum operation is done without concurency.";
CircTrappedIons::error = "`1`";
CircSiliconDelft::error = "`1`";
CircRydbergHub::error = "`1`";
Serialize::error = "`1`";

(* arrangement options *)
Parallel::usage = "Parallel options in arrangement of gates: False, Default, All. False: serial, default: parallel according to the device specification, and All: full quantum parallel";
MapQubits::usage = "Options in the CircTrappedIons[] that maps the local qubits (\[Rho]A) into the total qubits of the total (large) density matrix (\[Rho]AB).";

(** Custom gates, used in Aliases **)
BeginPackage["`CustomGates`"];

	SWAPLoc::usage = "Swap the spatial locations of two qubits. This is implemented in neutral atoms and trapped ions.";	
	SWAPLoc::error = "`1`";
	SWAPLoc::warning = "`1`";
	
	ShiftLoc::usage = "ShiftLoc[v] the physical coordinate of a qubit by some vector v. This is implemented in neutral atoms.";
	ShiftLoc::error = "`1`";
	ShiftLoc::warning = "`1`";
	
	Wait::usage = "Wait[\[CapitalDelta]t] gate, doing nothing/identity operation for duration \[CapitalDelta]t.";

	CZ::usage = "Controlled-Z operation.";
	CRx::usage = "Conditional Rx[\[Theta]] rotation on the nuclear 13C NV-center qubit, conditioned on the electron spin state.";
	CRy::usage = "Conditional Ry[\[Theta]] rotation on the nuclear 13C NV-center qubit, conditioned on the electron spin state.";
	CRot::usage = "Conditional rotation in Silicon spin qubit; this effectively implements CNOT gate.";
	Ent::usage = "Remote entanglement operation on multi-node trapped ions.";
	Splz::usage = "Splz[node, zone_destination]. Split a string of ions in a zone of a trapped-ion Oxford device";
	Shutl::usage = "Shutl[node,zone_destination]. Shuttle the qubit(s) to the destination zone";
	Comb::usage = "Comb[node, zone_destination]. Combine a string of ions to a zone of a trapped-ion Oxford device";
	PSW::usage = "PSW[\[Theta]], parameterised swaps";
	SRot::usage = "Single qubit gate in a driven Rydberg qubit via two-photon Raman transition. Usage: SRot[\[Phi],\[CapitalDelta],t] where \[Phi] is laser phase, \[CapitalDelta] is detuning, t is laser duration.";
	Init::usage = "Initialise qubit to its fiducial state; mostly to state \!\(\*TemplateBox[{\"0\"},\n\"Ket\"]\)";
	ZZ::usage = "The siZZle gate on a Superconducting device. Implemented by Exp[-(i\[Theta]/2) ZZ].";
	ZX::usage = "The cross-resonance gate on a Superconducting device. Implemented by Exp[-(i\[Theta]/2) ZX].";
	
	MPar::usage = "Perform parity measurement for two qubits that projects them into even (00,11) subspace and odd (01,10) subspace. In the case of Silicon qubit, state 01 decays to 10.";
	MPar::error="`1`";
	
EndPackage[]

(** All parameters used in virtual device configuration **)
BeginPackage["`ParameterDevices`"];

	AtomLocations::usage = "Three-dimensional physical locations of each atom/qubit.";
	Anharmonicity::usage = "The anharmonicity in the Superconducting device capturing the capacitor property.";
	BFProb::usage = "Probability of bit-flip error";
	BField::usage = "The electromagnetic field strength in the z-direction from the lab reference with unit Tesla.";
	BlockadeRadius::usage = "Short-range dipole-dipole interaction of Rydberg atoms in \[Mu]s. This allows multi-qubit gates.";
	Connectivity::usage = "Show the connectivity graph of a Superconducting qubit device, where the arrow show possible direction of the cross-resonant ZX gates.";
	DurTwoGate::usage = "Duration of two qubit gates with rotation of \[Pi]";
	DurMeas::usage = "Duration of measurement";
	DurInit::usage = "Duration of initialisation";
	DurRead::usage = "Readout duration in \[Mu]s";
	DurShuffle::usage = "Duration to shuffle location of ions";
	DurMove::usage = "Duration for physically moving operation in Trapped Ions such as Splz and Comb.";
	DurRxRy::usage = "Duration to run the rotation gates Rx and Ry, that is fixed regardless the angle.";
	DurZX::usage="Duration of the resonance ZX gate on the superconducting qubits that is fixed regardless the angle.";
	DurZZ::usage="Duration of the siZZle ZZ gate on the superconducting qubits that is fixed regardless the angle.";
	EFSingleXY::usage = "Error fraction/ratio, {depolarising, dephasing} of the single qubit X and Y rotations. Sum of the ratio must be 1 or 0 (off).";
	EFSingle::usage = "Error fraction/ratio, {depolarising, dephasing} of the single qubit X, Y, and Z rotations. Sum of the ratio must be 1 or 0 (off).";
	EFTwo::usage = "Error fraction/ratio, {depolarising, dephasing} of controlled-rotation. Sum of the ratio must be 1 or 0 (off).";
	EFCZ::usage = "Error fraction/ratio, {depolarising, dephasing} of the controlled-Z gates. Sum of the ratio must be 1 or 0 (off).";
	EFInit::usage = "Error fraction/ratio {depolarising, dephasing} of the initialisation gate.  Sum of the ratio must be 1 or 0 (off).";
	EFCRot::usage = "Error fraction/ratio, {depolarising, dephasing} of the controlled-Rx and -Ry gates. Sum of the ratio must be 1 or 0 (off).";
	EFEnt::usage = "Error fraction/ratio {depolarising, dephasing} of remote entanglement.  Sum of the ratio must be 1 or 0 (off).";
	EFRead::usage = "Error fraction/ratio of {depolarising,dephasing} of the readout. Sum of the ratio must be 1 or 0 (off)." ;
	ExchangeRotOn::usage = "Maximum interaction j on the passive qubit crosstalk when applying CZ gates; The noise form is C[Rz[j.\[Theta]]] It must be a square matrix with size (nqubit-2)x(nqubit-2).";
	ExchangeRotOff::usage = "Crosstalks error C-Rz[ex] on the passive qubits when not applying two-qubit gates.";
	ExcitedInit::usage = "The probability/fraction of the population excited in the thermal state. This is the same initialisation state.";
	ExchangeCoupling::usage = "The exchange coupling strength of resonators in Superconducting devices";
	FidCRot::usage = "Fidelity of conditional rotation in NV-center obtained by dynamical decoupling and RF pulse or the Controlled-X180 in the Silicon qubits that is used in readout/measurement.";
	FidSingleXY::usage = "Fidelity(ies) of single Rx[\[Theta]] and Ry[\[Theta]] rotations obtained by random benchmarking.";
	FidSingleZ::usage = "Fidelity(ies) of single Rz[\[Theta]] rotation obtained by random benchmarking.";
	FidSingle::usage = "Fidelity(ies) of single rotations: Rx[\[Theta]], Ry[\[Theta]], Rz[\[Theta]] obtained by random benchmarking.";
	FidTwo::usage = "Fidelity(ies) of two qubit gates obtained by random benchmarking.";
	FidEnt::usage = "Fidelity of remote entanglement operation.";
	FidMeas::usage = "Fidelity of measurement";
	FidInit::usage = "Fidelity of qubit initialisation";
	FidCZ::usage = "Fidelity(ies) of the CZ gates.";
	FreqSingleXY::usage = "Rabi frequency(ies) for the single X- and Y- rotations with unit MHz";
	FreqSingleZ::usage = "Rabi frequency(ies) for the single Z- rotations with unit MHz";
	FreqCZ::usage = "Rabi frequency(ies) for the CZ gate with unit MHz.";
	FreqEnt::usage = "Frequency of remote entanglement.";
	FreqCRot::usage = "Frequency of conditional rotation in NV-center obtained by dynamical decoupling and RF pulse.";
	FreqTwoGate::usage = "Resonant frequency of two qubit gates";
	FreqWeakZZ::usage = "Frequency of coherent cross-talk noise in form of ZZ-coupling that slowly entangle the qubits.";
	FidRead::usage = "Readout fidelity";
	GlobalField::usage = "Global magnetic field fluctuation in NV-center due to C13 bath. It causes dephasing on the electron in 4\[Mu]s.";
	HeatFactor::usage = "The constant >=1 that enhance dephasing process. This happends in neutral atoms when moving the atoms.";
	LossAtoms::usage = "Device key in the RydbergHub device that identifies atoms lost to the environment.";
	LossAtomsProbability::usage = "Device key in the RydbergHub that identifies probability of atoms lost to the environment due to repeated measurement.";
	(* future feature
	MoveSpeed : "The speed of moving atom \[Mu]m/\[Mu]s in Neutral Atom systems. This will affect the heat factor";
	*)
	Nodes::usage = "Entire nodes of a trapped ions system <|node1 -> number_of_qubits_1, ... |>";
	NIons::usage = "The total number of ions in a trapped ion device.";
	OffResonantRabi::usage = "Put the noise due to off-resonant Rabi oscillation when applying single qubit rotations.";
	QubitNum::usage = "The number of physical active qubits for computations.";
	QubitFreq::usage = "The fundamental qubit frequency for each qubit with unit MHz.";
	QMap::usage = "Show maps from nodes in trapped ions to the actual emulated qubits";
	ProbLeakInit::usage = "Leakage probability in the Rydberg initialisation. The noise is decribed with non-trace-preserving map.";
	ProbLeakCZ::usage = "Leakage probability in executi multi-controlled-Z.";
	ProbLossMeas::usage = "Probability of phyiscal atom loss due to measurement.";
	ProbLeakMove::usage = "Probability of leakage in the process of moving atoms.";
	ProbBFRot::usage = "Assymetric Bit-flip probability on single rotation operation. {01->p1, 10->p2}";
	RydbergRabiFreq::usage = "Rydberg Rabi frequency";
	RabiFreq::usage = "The Rabi frequency frequency in average or on each qubit with unit MHz.";
	ScatProb::usage = "Scattering probability";
	ShowNodes::usage = "Draw all Ions on every nodes within the zones";
	StdPassiveNoise::usage = "Set to True/False. Use the standard passive noise that involves T1, T2 or T2s inputs.";
	T1::usage = "T1 duration(s) in \[Mu]s. Exponential decay time for the state to be complete mixed.";
	T2::usage = "T2 duration(s) in \[Mu]s. Exponential decay time for the state to be classical with echo applied.";
	T2s::usage = "T2* duration(s) in \[Mu]s. Exponential decay time for the state to be classical.";
	UnitLattice::usage = "The unit lattice AtomLocations in \[Mu]s. This gives access to internal device parameter in RydbergHub device.";
	VacLifeTime::usage = "The lifetime of the qubit array is limited by its vacuum lifetime, where T1=VacTime/Nqubits.";
	ZZPassiveNoise::usage = "The swtich for ZZ interaction passive noise in the Superconducting device.";
	BField::error = "`1`";
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
	QubitNum::error="`1`";
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
	
	(*Validate expression, throw error if false*)
	validate[value_, expr_, err_, msg_, format_:Identity] :=
	    (
	        If[expr[format[value]],
	            format[value]
	            ,
	            Throw[Message[err::error, msg]]
	        ]
	    )
	(* check if it's an association with length len. *)
	checkAss[ass_, len_] :=
	    AssociationQ[ass] && Length[ass] === len && And @@ NumberQ /@ Values @ ass

	checkAss[ass_, len_, f_] :=
	    AssociationQ[ass] && Length[ass] === len && And @@ f /@ Values @ ass
	(* convert number to association *)
	num2Ass[arg_Real, len_Integer] :=
	    <|Table[i -> arg, {i, 0, -1 + len}]|>
	
	num2Ass[arg_Integer, len_Integer] :=
	    <|Table[i -> arg, {i, 0, -1 + len}]|>
	
	num2Ass[arg_Association, len_Integer] := 
		arg
	
	(* TODO *)
	numass[len_]:="not a number or association of numbers with length "<>ToString[len]
	fidass[len_]:="not a fidelity number or association of fidelities with length "<>ToString[len]
	
	(* convenient functions *)
	
	(* check if it is a probability number *)
	isProbability[p_ ]:=If[0 <= p <= 1, True, False]
	
	(*
		fid2DepolDeph[fidelity, {ratio.Depol,ratio.Deph}, nqubits, error_variable, is_average_fidelity:True]. 
		Return the parameter for depolarizing and dephasing noise that gives the total fidelity totFid [0,1], 
		where totFid is the average fidelity obtained from random benchmarking (set avgfid:False if it's the worst fidelity)
		correction is a constant adjusted to params of dephasing to get average fidelity from the worst fidelity.
	*)
	fid2DepolDeph[totfid_, errratio_, nqerr_, errval_, avgfid_:True] := Module[
		{sol, rate, pdepol, pdeph, entfid, d, p},
		(*set to perfect fidelity of zero error *)
		If[(totfid == 1 || 0 == Total @ errratio),
			{pdepol, pdeph} = {0, 0}
			,
			(* 1-qubit error case *)
			If[nqerr === 1,
				sol = If[avgfid,
					(* estimate parameters from entanglement fidelity *)
					entfid = (3 totfid - 1)/2; Solve[1 - rate * errratio[[1]] - rate * errratio[[2]] + 4/3 errratio[[1]] * errratio[[2]] * rate^2 == entfid,
						 {rate}]
					, 
					(*estimate parameters from the worst fidelity: from + state *)
					Solve[1 - 2 errratio[[1]] * rate/3 - errratio[[2]] * rate + 4 * errratio[[1]] * errratio[[2]] * rate^2/3 == totfid, {rate}]
				];
				          
				 (* check validity of numbers and set it to the worst parameter if exceeds *)
				{pdepol, pdeph} = errratio * Min @ Abs[rate /. sol];
				
				If[pdepol > 0.75, Message[errval::warning, StringForm["(warning) fidelity is too low; 1-qubit depolarization parameter is ``. Set it to 3/4.",
					 pdepol]]; pdepol = 0.75;];
				If[pdeph > 0.5, Message[errval::warning, StringForm["(warning) fidelity is too low; 1-qubit dephasing parameter is ``. Set it to 1/2.",
					 pdeph]]; pdeph = 0.5;];
				,
				(* 2-qubit error case*)
				sol = If[avgfid,
					(*estimate parameters from entanglement fidelity*)
					entfid = (5 * totfid - 1)/4; Solve[1 - rate * errratio[[1]] - rate * errratio[[2]] + 16/15 * rate^2 * errratio[[1]] * errratio[[2]] == entfid, {rate}]
					,
					(*estimate parameters from the worst fidelity: from ++ state   *)
					Solve[1 - errratio[[2]] * rate - errratio[[1]] * rate * 4/5 + errratio[[1]] * errratio[[2]] * rate^2 * 16/15 == totfid, {rate}]
				];
				
				(* check validity of numbers, set it to the worst if exceeds *)
				{pdepol, pdeph} = errratio * Min @ Abs[rate /. sol];
				If[pdepol > 15/16, Message[errval::warning, StringForm["(warning) Fidelity might be too low; 2-qubit depolarization parameter is ``. Set it to 15/16.",
					 pdepol]]; pdepol = 15/16;];
				If[pdeph > 3/4, Message[errval::warning, StringForm["(warning) Fidelity might be too low; 2-qubit dephasing parameter is ``. Set it to 3/4.",
					 pdeph]]; pdeph = 3/4;];
			];
		];
		{pdepol, pdeph}
	]

	(*
	Obtain the error parameters given entanglement fidelity. This works only for the trapped ion network fidelity.
	This is obtained by calculating <\[CapitalPsi]^+|Subscript[Depol, 0,1](d)Subscript[Deph, 0,1](p)|\!\(
\*SuperscriptBox[\(\[CapitalPsi]\), \(+\)]
\*SubscriptBox[\(>\), \(0, 1\)]\)
	*)
	entFid2DepolDeph[entfid_, errratio_, errval_:FidEnt] := Module[
		{pdepol, pdeph, sol, rate}
		,
		sol = Solve[1 - rate * errratio[[1]] * 4/5 - rate * errratio[[2]] * 2/3 + rate^2 * errratio[[1]] * errratio[[2]] * 32/45 == entfid,
			 {rate}];
		{pdepol, pdeph} = errratio * Min @ Abs[rate /. sol];
		If[pdepol > 15/16,
			Message[errval::warning, StringForm["(warning) Fidelity might be too low; 2-qubit depolarization parameter is ``. Set it to 15/16.",
				 pdepol]]; pdepol = 15/16;
		];
		If[pdeph > 3 / 4,
			Message[errval::warning, StringForm["(warning) Fidelity might be too low; 2-qubit dephasing parameter is ``. Set it to 3/4.",
				 pdeph]]; pdeph = 3/4;
		];
		{pdepol, pdeph}
	]
(** Error channels expressed as a list of matrices **)

	(* 1-qubit bit-flip error *)
	bitFlip1[fid_] :=
		With[{e = 1 - fid},
			{Sqrt[1 - e] * {{1, 0}, {0, 1}}, Sqrt[e] * {{0, 1}, {1, 0}}}
		]
	
	(* 2-qubit bit-flip error *)
	bitFlip2[fid_] :=
		With[{e = 1 - fid},
			{Sqrt[1 - e] * {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}, 
			Sqrt[e/3] * {{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}}, 
			Sqrt[e/3] * {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}, {0, 1, 0, 0}}, 
			Sqrt[e/3] * {{0, 0, 0, 1}, {0, 0, 1, 0}, {0, 1, 0, 0}, {1, 0, 0, 0}}}
		]

	(* 
	Generalised amplitude damping (Nielsen&Chuang, P.382).
	Common description of T1 decay, by setting gamma(t)=1-Exp[-t/T1],
	which describes the shrinking Bloch sphere into the ground state zero,
	with probability p.
	*)
	(* Rho_infty = p|0X0|+(1-p)|1X1|; this returns the operator *)
	Subscript[gAmp, q_][gamma_, p_] :=
		Subscript[Kraus, q][
			{
				Sqrt[p]*{{1, 0}, {0, Sqrt[1 - gamma]}}
				,
				Sqrt[p]*{{0, Sqrt[gamma]}, {0, 0}}
				,
				Sequence @@
					If[p < 1,
						{Sqrt[1 - p]*{{Sqrt[1 - gamma], 0}, {0, 1}}, Sqrt[1 - p] * {{0, 0}, {Sqrt[gamma], 0}}}
						,
						{}
					]
			}
		]
		
	(* Off-resonant Rabi oscillation. 
		omega: Rabi frequency
		delta: detuning frequency
		dt: time duration of the pulse
		*)
	offresonantRabi[omega_, delta_, dt_] :=
		With[{omegaR = Sqrt[omega^2 + delta^2]},
			E^(I * delta * dt/2)*
					{{Cos[omegaR * dt/2] - I * Sin[omegaR * dt/2] * delta/omegaR, -I * Sin[omegaR * dt/2] * omega/omegaR},
	                 {-I * Sin[omegaR * dt/2]*omega/omegaR, Cos[omegaR * T/2] + I * Sin[omegaR * T/2] * delta/omegaR}}
		]
(** Custom gates definitions **)
	(* Unitary driven by laser detuning *)												
	unitaryDetuning[phi_, delta_, dt_, rabi_] :=
		Module[{vdelta, vt, omega},
			FullSimplify[
				{{Cos[omega * vt/2] - I * Sin[omega * vt/2] * vdelta/omega, -I * E^(I*phi) Sin[omega * vt/2]* rabi/omega}, 
				{-I * E ^ (-I * phi) * Sin[omega * vt/2] * rabi/omega, Cos[omega vt/2] + I * Sin[omega * vt/2] * delta/omega}} //. 
				{vdelta -> delta, omega -> Sqrt[rabi^2 + vdelta^2], vt -> dt}
				 , {rabi >= 0, vdelta >= 0}]
		]
																																															
	(************  VIRTUAL_DEVICES  *************)		
																																														
(* DEVICE_TOY *)
	ToyDevice[OptionsPattern[]] := With[
	{
		qubitnum = OptionValue @ QubitNum, 
		t1 = OptionValue @ T1, 
		t2 = OptionValue @ T2,
		stdpassivenoise = OptionValue @ StdPassiveNoise, 
		fidsingle = OptionValue @ FidSingle, 
		fidtwo = OptionValue @ FidTwo, 
		efsingle = OptionValue @ EFSingle, 
		eftwo = OptionValue @ EFTwo, 
		rabifreq = OptionValue @ RabiFreq, 
		twogatefreq = OptionValue @ FreqTwoGate
	},
		
			Module[{deltaT, erone, ertwo},
				erone = fid2DepolDeph[fidsingle, efsingle, 1, FidSingle, True];
				ertwo = fid2DepolDeph[fidtwo, eftwo, 2, FidTwo, True];
				<|
					(*no hidden qubits/ancilla here *)
					DeviceType -> "Toy"
					,
					OptionsUsed -> {
						QubitNum -> qubitnum,
						T1 -> t1,
						T2 -> t2,
						StdPassiveNoise -> stdpassivenoise, 
						FidSingle -> fidsingle, 
						FidTwo -> fidtwo, 
						EFSingle -> efsingle, 
						EFTwo -> eftwo, 
						RabiFreq -> rabifreq, 
						FreqTwoGate -> twogatefreq
					}
					,
					DeviceDescription -> "Toy device with " <> ToString[qubitnum] <> "-qubits arranged as a linear array with nearest-neighbor connectivity."	
					,
					NumAccessibleQubits -> qubitnum
					,
					NumTotalQubits -> qubitnum
					,
					Aliases -> {
						(* parameterised swap *)
						Subscript[PSW, p_, q_][theta_] :> Subscript[U, p, q][
							{{1, 0, 0, 0}, {0, E^((I theta)/2) Cos[theta/2], -I E^((I theta)/2) Sin[theta/2], 0},
							 {0, -I E ^ ((I theta)/2) Sin[theta/2], E ^ ((I theta)/2) Cos[theta/2], 0}, {0, 0, 0, 1}}]
						 }
					,
					Gates ->
						{
							(* Singles *)
							
							Subscript[Rx, q_][theta_] :>
								<|
								NoisyForm -> Flatten @ {Subscript[Rx, q][theta], Subscript[Depol, q][erone[[1]]], Subscript[Deph, q][erone[[2]]]},
								GateDuration -> Abs[theta]/rabifreq
								|>
							,
							Subscript[Ry, q_][theta_] :> 
								<|
								NoisyForm -> Flatten @ {Subscript[Ry, q][theta], Subscript[Depol, q][erone[[1]]], Subscript[Deph, q][erone[[2]]]},
								GateDuration -> Abs[theta]/rabifreq
								|>
							,
							Subscript[Rz, q_][theta_] :> 
								<|
								NoisyForm -> Flatten @ {Subscript[Rz, q][theta], Subscript[Depol, q][erone[[1]]], Subscript[Deph, q][erone[[2]]]}, 
								GateDuration -> Abs[theta]/rabifreq
								|>
							,
							(* Twos *)
							Subscript[C, p_][Subscript[Rx, q_][theta_]] /; Abs[q - p] === 1 :> 
								<|
								NoisyForm -> {Subscript[C, p][Subscript[Rx, q][theta]], Subscript[Depol, p, q][Min[ertwo[[1]] * Abs[theta/\[Pi]], 15/16]], Subscript[Deph, p, q][Min[ertwo[[2]] * Abs[theta/\[Pi]], 3/4]]}, 
								GateDuration -> Abs[theta]/twogatefreq
								|>
							,
							Subscript[C, p_][Subscript[Ry, q_][theta_]] /; Abs[q - p] === 1 :> 
								<|
								NoisyForm -> {Subscript[C, p][Subscript[Ry, q][theta]], Subscript[Depol, p, q][Min[ertwo[[1]] * Abs[theta/\[Pi]], 15/16]], Subscript[Deph, p, q][Min[ertwo[[2]] * Abs[theta/\[Pi]], 3/4]]}, 
								GateDuration -> Abs[theta]/twogatefreq 
								|>
							,
							Subscript[C, p_][Subscript[Rz, q_][theta_]] /; Abs[q - p] === 1 :> 
								<|
								NoisyForm -> {Subscript[C, p][Subscript[Rz, q][theta]], Subscript[Depol, p, q][Min[ertwo[[1]] * Abs[theta/\[Pi]], 15/16]], Subscript[Deph, p, q][Min[ertwo[[2]] * Abs[theta/\[Pi]], 3/4]]}, 
								GateDuration -> Abs[theta]/twogatefreq
								|>
							, 
							(* parameterised swap*)
							Subscript[PSW, p_, q_][theta_] /; Abs[q - p] === 1 :> 
								<|
								NoisyForm -> {Subscript[PSW, p, q][theta], Subscript[Depol, p, q][ertwo[[1]]], Subscript[Deph, p, q][ertwo[[2]]]}, 
								GateDuration -> Abs[theta]/twogatefreq
								|>
						}
					,
					(* Declare that deltaT will refer to the duration of the current gate/channel. *)
					DurationSymbol -> deltaT
					,
					(* Passive noise *)
					Qubits :>
						{
							q_Integer :>
								<|
									PassiveNoise ->
										(* standard passive noise on *)
										If[stdpassivenoise,
											{Subscript[Depol, q][0.75 (1 - E ^ (-deltaT / t1))], Subscript[Deph, q][0.5 (1 - E ^ (-deltaT / t2))]}
											,
											{}
										]
								|>
						}
				|>
			]
		]

(* DEVICE_SUPERCONDUCTING_QUBITS *)

	(* 
	construct graph that shows connectivity of the corss-resonance gate (ZX) in SQC based on the value of qubit frequencies
	 *)
	graphConnectivityFromQubitFreq[nqubit_, coupling_, qubitFreq_] :=
	Module[
		{g, dedge}
		,
		(*directed edges, control->target*)
		dedge =
			If[qubitFreq[#[[1]]] > qubitFreq[#[[2]]],
					#[[1]] -> #[[2]]
					,
					#[[2]] -> #[[1]]
				]& /@ Keys[coupling];
		(* basic graph *)
		g = Graph[Range[0, nqubit - 1], dedge, VertexWeight -> Values @ qubitFreq, EdgeWeight -> Values @ coupling];
		
		(* the final graph with a visualisation format*)
		DirectedGraph[g, VertexSize -> 0.5, BaseStyle -> {15, Bold, FontFamily -> "Serif"}, 
			VertexLabels -> {v_ :> Placed[{AnnotationValue[{g, v}, VertexWeight], "Q" <> ToString[v]}, {Center, {After, Above}}]},
			EdgeLabels -> "EdgeWeight", GraphLayout -> "SpringEmbedding", EdgeStyle -> {Thick}, VertexStyle -> {Yellow, EdgeForm[None]}
		 ]
	]
	
	(* 
	device description of superconducting qubits 
	*)
	SuperconductingHub[OptionsPattern[]] :=
	With[{
		qubitsnum = OptionValue @ QubitNum,
		t1 = OptionValue @ T1,
		t2 = OptionValue @ T2,
		excitedinit = OptionValue @ ExcitedInit,
		qubitfreq = OptionValue @ QubitFreq,
		exchangecoupling = OptionValue @ ExchangeCoupling,
		anharmonicity = OptionValue @ Anharmonicity,
		fidread = OptionValue @ FidRead,
		durmeas = OptionValue @ DurMeas,
		durrxry = OptionValue @ DurRxRy,
		durzx = OptionValue @ DurZX,
		durzz = OptionValue @ DurZZ,
		stdpassivenoise=OptionValue @ StdPassiveNoise,
		zzpassivenoise=OptionValue @ ZZPassiveNoise
		},
		
		(* some assertions to check the inputted parameters *)
		Catch @ If[CountDistinct[Values@qubitfreq]==Length@qubitfreq, 
			Throw @ Message[QubitFreq::error, "Every qubit frequency value given in QubitFreq must be distinct."]
		];
		
		Module[{ccv, ug, stdpn, zzInteraction, zzPN, lessNeighbor, zzon, passivenoise, deltaT, activeq, counter=0, init=True}
			,		
			(* cross-resonance connectivity *)
			ccv = graphConnectivityFromQubitFreq[qubitsnum, exchangecoupling, qubitfreq];
			
			(* undirected graph *)	
			ug = UndirectedGraph[ccv];
			
			(* track the active qubit for passive noice application *)
			activeq = <|Table[q -> False, {q, Range[0, qubitsnum - 1]}]|>;
			
			(* track if a there is ZZ-gate (siZZler) applied *)
			zzon = False;
			
			(* Free-induction T1- T2- passive noise decays *)	
			Subscript[stdpn, q_][dur_] := If[stdpassivenoise && \[Not]activeq[q] && dur > 0,
				{Subscript[gAmp, q][(1 - E ^ (-dur/t1[q])), 1 - excitedinit[q]], Subscript[Deph, q][0.5 (1 - E ^ (-dur / t2[q]))]}
				,
				{}
			];
			
			(* Fixed ZZ-interaction on passive noise, where c is the first/control and t is the later/target *)		
			Subscript[zzInteraction, c_, t_] := With[{
				deltact = Abs[qubitfreq[c] - qubitfreq[t]]
				,
				alphat = anharmonicity[t]
				,
				alphac = anharmonicity[c]
				,
				jc = If[KeyExistsQ[exchangecoupling, c \[UndirectedEdge] t],
					exchangecoupling[c \[UndirectedEdge] t]
					,
					exchangecoupling[t \[UndirectedEdge] c]
				]
			},
				R[jc^2 (1 / (deltact - alphat) - 1 / (deltact + alphat)), Subscript[Z, c] Subscript[Z, t]]
			];
			
			(* list neighbor qubits with ordering less than q in graph g. *)		
			lessNeighbor[q_] := List @@@ Select[EdgeList[NeighborhoodGraph[ccv, q]], q > First @ DeleteElements[List @@ #, {q}]&];
			
			(* zz-crosstalk as passive noise which is off when a siZZle is on *)
			Subscript[zzPN, q_][dur_] := Module[{ng = lessNeighbor[q], noise},
				noise = If[zzpassivenoise,
					(* a siZZle gate is active in a parallel column *)
					If[\[Not]zzon,
						(
							If[dur > 0 && And[\[Not]activeq[#[[1]]], \[Not]activeq[#[[2]]]],
								Subscript[zzInteraction, #[[1]], #[[2]]]
								,
								{}
							]
						)& /@ ng
						,
						{}
					]
					,
					{}
				];
				noise
			];
					
		<|
			(* store the option that is initially used here *)
			OptionsUsed -> 
			{
				QubitNum -> qubitsnum,
				T1 -> t1,
				T2 -> t2,
				ExcitedInit -> excitedinit,
				QubitFreq -> qubitfreq,
				ExchangeCoupling -> exchangecoupling,
				Anharmonicity -> anharmonicity,
				FidRead -> fidread,
				DurMeas -> durmeas,
				DurRxRy -> durrxry,
				DurZX -> durzx,
				DurZZ -> durzz,
				StdPassiveNoise -> stdpassivenoise,
				ZZPassiveNoise -> zzpassivenoise
			}
			,
			DeviceType->"Superconducting"
			,
			DeviceDescription -> ToString[qubitsnum]<>"-qubit of Superconducting transmon qubits based on Josephson junctions"
			,
			NumAccessibleQubits -> qubitsnum
			,
			NumTotalQubits -> qubitsnum
			,
			Connectivity -> ccv
			,
			Aliases -> 
			{		
				Subscript[ZZ, p_,q_] :> R[\[Pi]/2,Subscript[Z, p] Subscript[Z, q]]
				,
				Subscript[ZX, p_,q_] :> R[\[Pi]/2,Subscript[Z, p] Subscript[X, q]]
				,
				(* in practice, it is a decay to the thermal state *)
				Subscript[Init, q_] :> Nothing
				,
				Subscript[Wait, q_] :> Nothing
			},	
			Gates ->{
				(* apply generalised amplitude damping to describe the decay to the thermal state which works only in the beginning *)
				Subscript[Init, q_] /; (init == True) :> 
				<|
					NoisyForm -> {Subscript[Init, q], Subscript[gAmp, q][1, 1 - excitedinit[q]]},
					GateDuration -> 0,
					UpdateVariables -> Function[
							activeq[q] = True;
						]
				|>
				,
				(* TODO: assertion for the readout at the end only *)
				Subscript[M, q_] :>
				<|
					(* depolarise up the final result as well *)
					NoisyForm -> {Subscript[gAmp, q][(1 - E ^ (-durmeas/t1[q]))],Subscript[Depol, q][1 - fidread[q]], Subscript[M, q], Subscript[Depol, q][1-fidread[q]]},
					GateDuration -> durmeas,
					UpdateVariables -> Function[
							activeq[q] = True;
							init = False;
						]
				|>
				,
				(* doing nothing, is equivalent to being passive *)
				Subscript[Wait, q_][dur_]:>
				<|
					NoisyForm -> Flatten @ {Subscript[stdpn, q][dur], Subscript[zzPN, q][dur]},
					GateDuration -> dur,
					UpdateVariables -> Function[ activeq[q] = False ]
				|>
				,
				
				(* single-qubit gates  *)
					Subscript[Rx, q_][theta_] /; If[NumberQ@theta, And[(-Pi <= theta <= Pi), Abs[theta] > 0], True] :>
					<|
						NoisyForm -> {Subscript[Rx, q][theta]},
						GateDuration -> durrxry,
						UpdateVariables -> Function[
							activeq[q] = True;
							init = False;
						]
					|>
					,
					Subscript[Ry, q_][theta_] /; If[NumberQ@theta, And[(-Pi <= theta <= Pi), Abs[theta] > 0], True] :>
					<|
						NoisyForm -> {Subscript[Ry, q][theta]},
						GateDuration -> durrxry,
						UpdateVariables -> Function[
							activeq[q] = True;
							init = False;
						]
					|>
					,	
					(* virtual gate *)			
					Subscript[Rz,q_][theta_] :>
					<|
						NoisyForm -> Flatten@{Subscript[Rz, q][theta]},
						GateDuration -> 0
					|>
					,				
					(** two-qubit gates **)
					(* siZZle gate *)
					Subscript[ZZ, p_, q_] /; EdgeQ[ug, p\[UndirectedEdge]q] :>
					<|
						NoisyForm -> {Subscript[ZZ, p, q]}, 
						GateDuration -> durzz,
						UpdateVariables -> Function[
							init = False;
							zzon = True;	
							activeq[q] = True;
							activeq[p] = True;
							]
					|>
					,
					(* cross-resonance gate *)
					Subscript[ZX, p_, q_] /; EdgeQ[ccv, p\[DirectedEdge]q] :>
					<|
						NoisyForm -> {Subscript[ZX, p, q]}, 
						GateDuration -> durzx,
						UpdateVariables -> Function[
							init = False;
							activeq[q] = True;
							activeq[p] = True;
							]
					|>		
			},
			
			(*TODO: option to activate/deactivate this. Re-initialized when invoking InsertCircuitNoise *)
			InitVariables -> 
				Function[
					activeq = Association[ # -> False & /@ Range[0, qubitsnum - 1]];
					zzon = False;
					init = True;
				]
			,
			
			(* Passive noise *)
			(* Declare that deltaT will refer to the duration of the current gate/channel. *)
			DurationSymbol -> deltaT
			, 
			Qubits :> {
				q_ :> <|
						(* reset the zzon here *)
						PassiveNoise -> Flatten @ {Subscript[stdpn, q][deltaT], Subscript[zzPN, q][deltaT]}
						,
						UpdateVariables -> Function[
							counter++; 
							If[counter == qubitsnum,	
								zzon = False;
								Table[activeq[j] = False, {j, Keys @ activeq}];
								counter = 0;
							]
						]	
						|>	
						
					}
			|>
		]
	]
	
	
(* DEVICE_NVCENTER_DELFT *)
	NVCenterDelft[OptionsPattern[]] := With[
	{
		qubitnum = OptionValue @ QubitNum,
		t1 = OptionValue @ T1,
		t2 = OptionValue @ T2,
		freqcrot = OptionValue @ FreqCRot,
		freqsinglexy = OptionValue @ FreqSingleXY,
		freqsinglez = OptionValue @ FreqSingleZ,
		fidcrot = OptionValue @ FidCRot,
		fidsinglexy = OptionValue @ FidSingleXY,
		fidsinglez = OptionValue @ FidSingleZ,
		efsinglexy = OptionValue @ EFSingleXY,
		efcrot = OptionValue @ EFCRot,
		fidinit = OptionValue @ FidInit,
		fidmeas = OptionValue @ FidMeas,
		durmeas = OptionValue @ DurMeas,
		durinit = OptionValue @ DurInit,
		freqweakzz = OptionValue @ FreqWeakZZ
	},
	
	Module[
		{dt, stdpn, passivenoise, ersinglexy, ersinglez, ercrot, nuclearq, weakzz}
		,
		(* standard passive noise *) 
		stdpn[q_, dur_] := 
			{Subscript[Depol, q][0.75 (1 - E^(-dur/t1[q]))], Subscript[Deph, q][0.5 (1 - E^(-dur/t2[q]))]};
		(* 
		Note: cross-talk ZZ-coupling in order of Hz on passive noise, Phase entanglement among nuclear spins 
		Subscript[C, n1][Subscript[Rz, n2][dt]], weak rotation for all combinations of n1,n2, few Hz.
		implement: Exp[-i dur ZZ]
		*)	
		weakzz[q_, dur_] := 
			R[dur * Pi/freqweakzz, Subscript[Z, #] Subscript[Z, #2]]& @@@ Subsets[DeleteElements[Range[1, qubitnum - 1], {q}], {2}];
	
		passivenoise[q_, dur_] :=
			If[NumberQ @ freqweakzz,
				Flatten @ {stdpn[q, dur], weakzz[q, dur]}
				,
				stdpn[q, dur]
			];
	
		(* error parameters (depolarising and dephasing) *)
		ersinglexy = fid2DepolDeph[#, efsinglexy, 1, FidSingleXY]& /@ fidsinglexy;	
		
		ersinglez = fid2DepolDeph[#, {0, 1}, 1, FidSingleZ]& /@ fidsinglez;	
		
		ercrot = fid2DepolDeph[#, efcrot, 2, FidCRot]& /@ fidcrot;
	
	<|
		DeviceDescription ->"One node of an NV center, where qubit 0 is the electronic spin. It has start connectivity with qubit 0 at the center.",
		
		(* The number of accessible qubits. This informs the qubits that a user's circuit can target *)
		NumAccessibleQubits -> qubitnum,	
		
		(* The total number of qubits which would be needed to simulate a circuit prescribed by this device.
		 * This is > NumAccessibleQubits only when the spec uses hidden qubits for advanced noise modelling *)
		NumTotalQubits -> qubitnum,
	
		(* Aliases are useful for declaring custom events. At this stage the specification is noise-free (but see later) *)
	Aliases -> {
		Subscript[Init, q_] :>  Nothing
		,
		Subscript[Wait, q__][t_] :> Nothing
		,
		Subscript[CRx, e_, n_][theta_] :> Subscript[U, e, n][
				{{Cos[theta/2], 0, -I Sin[theta/2], 0}, {0, Cos[theta/2], 0, I Sin[theta/2]}, 
				{-I Sin[theta/2], 0, Cos[theta/2], 0}, {0, I Sin[theta/2], 0, Cos[theta/2]}}
			],
		Subscript[CRy, e_, n_][theta_] :> Subscript[U, e, n][
				{{Cos[theta/2], 0, -Sin[theta/2], 0}, {0, Cos[theta/2], 0, Sin[theta/2]},
				{Sin[theta/2], 0, Cos[theta/2], 0},{0, -Sin[theta/2], 0, Cos[theta/2]}}
			]
		}
		,	
	Gates ->
	{
	(* exclusively on ELECTRON SPIN, q===0 *)
			Subscript[Init, q_] /; q === 0  :>
			<|
				NoisyForm -> {Subscript[Damp, q][fidinit]}, 
				GateDuration -> durinit
			|>
			,
			Subscript[M, q_] /; q === 0 :> 
			<|
				NoisyForm -> {Subscript[X, 0], Subscript[Damp, 0][1-fidmeas], Subscript[X, 0], Subscript[M, 0]},
				GateDuration -> durmeas
			|>
			,
			(* Electron and nuclear spins *)
			(* A simple 'wait' instruction, useful for padding a circuit. The content inside the table should be the same as the passive noise section below (copy / paste it) *)
			Subscript[Wait, qubits__][dur_] :>
			<|
				NoisyForm -> stdpn[#, dur]& /@ Flatten[{qubits}],
				GateDuration -> dur
			|>
			,			
			(* Z-rotations are unconditional for both electron and nuclear spins *)
			Subscript[Rz, q_][theta_] :> 
			<|
				NoisyForm -> {Subscript[Rz, q][theta],Subscript[Deph, q][Min[ersinglez[q][[2]] * Abs[theta]/Pi,0.5]]},
				GateDuration -> Abs[theta]/freqsinglez[q]
			|>
			, 
			Subscript[Rx, q_][theta_] /; q === 0 :> 
			<|
				NoisyForm -> {Subscript[Rx, q][theta], Subscript[Depol, q][Min[ersinglexy[q][[1]] Abs[theta]/Pi, 0.75]], Subscript[Deph, q][Min[ersinglexy[q][[2]] Abs[theta]/\[Pi], 0.5]]},
				GateDuration -> Abs[theta]/freqsinglexy[q]
			|>
			,
			Subscript[Ry, q_][theta_]/; q === 0 :> <|
				NoisyForm -> {Subscript[Ry, q][theta], Subscript[Depol, q][Min[ersinglexy[q][[1]]Abs[theta]/Pi, 0.75]], Subscript[Deph, q][Min[ersinglexy[q][[2]] Abs[theta]/Pi, 0.5]]},
				GateDuration -> Abs[theta]/freqsinglexy[q]
			|>
			,
			(* NUCLEAR SPIN ROTATIONS CONDITIONED ON ELETRON SPIN. electron needs to be set at ms=-1: quite bad *)
			Subscript[Rx, q_][theta_] /; q > 0 :> 
			<|
				NoisyForm -> {Subscript[Rx, q][theta], Subscript[Depol, q][Min[ersinglexy[q][[1]]Abs[theta]/Pi, 0.75]], Subscript[Deph, q][Min[ersinglexy[q][[2]] Abs[theta]/Pi, 0.5]]},
				GateDuration -> Abs[theta]/freqsinglexy[q]
			|>
			,
			Subscript[Ry, q_][theta_] /; q > 0 :>
			<|
				NoisyForm -> {Subscript[Ry, q][theta], Subscript[Depol, q][Min[ersinglexy[q][[1]] Abs[theta]/Pi, 0.75]], Subscript[Deph, q][Min[ersinglexy[q][[2]] Abs[theta]/Pi, 0.5]]},
				GateDuration -> Abs[theta]/(freqsinglexy[q])
			|>
			,
			(* conditional rotations *)
			Subscript[CRx, e_, n_][theta_] /; (e == 0 && n > 0) :> 
			<|
				(* its noisy form depolarising the control and target qubits *)
				NoisyForm -> {Subscript[CRx, e, n][theta], Subscript[Depol,e, n][Min[ercrot[n][[1]] Abs[theta]/Pi, 15/16]], Subscript[Deph, e, n][Min[ercrot[n][[2]]Abs[theta]/Pi, 3/4]]},
				GateDuration -> Abs[theta]/(freqcrot[n]) 
			|>
			,
			Subscript[CRy, e_, n_][theta_] /; (e == 0 && n > 0):> 
			<|
				(* its noisy form depolarises the control and target qubits *)
				NoisyForm -> {Subscript[CRy, e, n][theta], Subscript[Depol, e, n][Min[ercrot[n][[1]] Abs[theta]/Pi, 15/16]], Subscript[Deph, e, n][Min[ercrot[n][[2]]Abs[theta]/Pi, 3/4]]},
				GateDuration -> Abs[theta]/(freqcrot[n])    
			|>		
		}	
		,
		(* Declare that dt will refer to the duration of the current gate/channel. *)
		DurationSymbol -> dt, 
	
		(* passive noise *)
		(* 
		Note 'globalField' is the unwanted (positive or negative) field offset for the present circuit; 
		this should be set on a per-run basis as in examples below. Could be augmented with e.g. a drifting function of t 
		other potential passive noise: Global magnetic field fluctuation due to ^13C bath 
		causes dephasing on the electron in 4 \[Mu]s
		external field is effectively perfect
		*)
			Qubits -> {
				q_ :> 
					<|
						PassiveNoise -> passivenoise[q, dt]
					|>	
			}
		|>
	
		]
	]
	
	
	(*  DEVICE_RYDBERGHUB  *)

	(* 
	Asymmetric bit-flip error where where one can assign the bit-flip of 0->1 and 1->0 separately. 
	Realised with amplitude damping and symmetric bitflip
	*)	
	Subscript[asymBitFlip, q_][p01_, p10_] := Module[
		{pbf = Min[p01, p10], pmax = Max[p01, p10], edamp, x}
		,
		 (*Here, we assume damping and bitflip are independent events. 
		P(damp or bf)=P(damp)+P(bf)-P(damp)P(bf)*)
		edamp = First[x /. Solve[x + pbf - x * pbf == pmax, {x}]];
		Which[(* more states are flipped to 1*)
			p01 > p10,
				Sequence @@ {Subscript[X, q], Subscript[Damp, q][edamp], Subscript[X, q], Subscript[Kraus, q][bitFlip1[1 - pbf]]}
			,
			p01 < p10,
				(*more states are flipped to 0*)Sequence @@ {Subscript[Damp, q][edamp], Subscript[Kraus, q][bitFlip1[1 - pbf]]}
			,
			True,
				(*equals*)Subscript[Kraus, q][bitFlip1[1 - pbf]]
		]
	]			

	(* 
	Plot for the neutral atoms configuration 
	*)
	SetAttributes[PlotAtoms, HoldFirst]
	PlotAtoms[rydbergdev_, opt : OptionsPattern[{PlotAtoms, Graphics}]] := With[{
		qulocs = rydbergdev @ AtomLocations, 
		lossqulocs = rydbergdev @ LossAtoms,
		blrad = rydbergdev @ BlockadeRadius, 
		unit = rydbergdev @ UnitLattice, 
		blockade = OptionValue @ ShowBlockade,
		showloss = OptionValue @ ShowLossAtoms,
		style = {ImageSize -> Medium, Frame -> True, Axes -> True}	
		},
		Which[
			2 === Length @ First @ Values @ qulocs,
			(*2D lattice*)
				Show[
					Sequence @@ Table[Graphics[{Cyan, Opacity[0.15], EdgeForm[Directive[Dashed, Orange]], Disk[unit * qulocs[b], blrad]}], {b, blockade}],
					Sequence @@ Table[Graphics[{Red, Disk[unit * v, 0.15]}], {v, Values @ qulocs}],
					Sequence @@ Table[Graphics[Text[k, ({0.1, 0.1} + qulocs[k]) * unit]], {k, Keys @ qulocs}],
					
					(* show the atoms lost to the environment at the last position *)
					If[showloss,
						Sequence @@ Table[Graphics[{Gray, Disk[unit * v, 0.15]}], {v, Values @ lossqulocs}]
						,
						Sequence @@ {}
					]
					,
					If[showloss,
						Sequence @@ Table[Graphics[Text[k, unit * ({0.15, 0.15} + qulocs[k])]], {k, Keys @ lossqulocs}]
						,
						Sequence @@ {}
					]
					,
					Evaluate @ FilterRules[{opt}, Options[Graphics]],
					Sequence @@ style,
					AxesLabel -> {"x", "y"}
				]
			,
			3 === Length @ First @ Values @ qulocs,
				Show[
					Sequence @@ Flatten @ Table[{Graphics3D[{Text[Style[k, Bold, White], unit * qulocs[k]]}], Graphics3D[{Red, Sphere[qulocs[k] * unit, 0.15]}]}, {k, Keys @ qulocs}],
					Sequence @@ Table[Graphics3D[{Cyan, Opacity[0.15], Sphere[unit * qulocs[b], blrad]}], {b, blockade}],
					(* show the atoms lost to the environment: the last position *)
					If[showloss,
						Sequence @@ Table[Graphics3D[{GrayLevel[0.5], Sphere[unit * v, 0.15]}], {v, Values @ lossqulocs}]
						,
						Sequence @@ {}
					]
					,
					If[showloss,
						Sequence @@ Table[Graphics3D[Text[Style[k, Bold, White], unit * qulocs[k]]], {k, Keys @ lossqulocs}]
						,
						Sequence @@ {}
					],
					Evaluate @ FilterRules[{opt}, Options[Graphics]],
					Sequence @@ style,
					AxesLabel -> {"x", "y", "z"}
				]
			,
			True
			,
			Message[PlotAtoms::error, "Only accepts 2D- or 3D-coordinates."];
			Return @ $Failed
		]
	]
		
		
	
	(* 
	The main virtual Rydberg function 
	*)																																																																		
	RydbergHub[OptionsPattern[]] := With[
		{
		qubitsnum = OptionValue @ QubitNum,
		atomlocations = OptionValue @ AtomLocations,
		t2 = OptionValue @ T2,
		vaclifetime = OptionValue @ VacLifeTime,
		rabifreq = OptionValue @ RabiFreq,
		unitlattice = OptionValue @ UnitLattice,
		blockaderad = OptionValue @ BlockadeRadius,
		probleakinit = OptionValue @ ProbLeakInit,
		probbfrot01 = OptionValue[ProbBFRot][01],
		probbfrot10 = OptionValue[ProbBFRot][10],
		durinit = OptionValue @ DurInit,
		fidmeas = OptionValue @ FidMeas,
		durmeas = OptionValue @ DurMeas,
		durmove = OptionValue @ DurMove,
		heatfactor = OptionValue @ HeatFactor,
		problossmeas = OptionValue @ ProbLossMeas,
		probleakcz = OptionValue @ ProbLeakCZ,
		qubits = Keys @ OptionValue @ AtomLocations
		},
	
		(* assertions *)
		If[Length @ atomlocations != qubitsnum,
				Message[QubitNum::error, "Missing or extra qubits in AtomLocations"]; Return @ $Failed];
		If[\[Not]isProbability[probleakinit],
				Message[ProbLeakInit::error, "Needs value within [0,1]"]; Return @ $Failed];
		If[\[Not]isProbability[problossmeas],
				Message[ProbLossMeas::error, "Needs value within [0,1]"]; Return @ $Failed];
		If[\[Not]isProbability[probleakcz],
				Message[ProbLeakCZ::error, "Needs value within [0,1]"]; Return @ $Failed];
	
		Module[
		{\[CapitalDelta]t,  lossatomprob, globaltime, stdpn, movenoise, t1, atomlocs, lossatomlocs, distloc, blockadecheck, circorloss, legshift}
		,
		(* record the location of each atom *)
		atomlocs = atomlocations;
		
		(* record the last position of the loss atom; this records the atoms that are lost to the environment at the same time *)
		lossatomlocs = <||>;
		
		(* Track loss probability of each atom. After measurement, a dice is thrown to decide if atom is still there *)
		lossatomprob = <|Table[k -> 0, {k, Keys @ atomlocations}]|>;
		
		(* Check if the atoms are there, otherwise send a warning and just apply complete depolarising noise *)
		circorloss[gate_, circ_, q__]:= With[{lostq = Intersection[Keys @ lossatomlocs, {q}]},
				If[Length @ lostq > 0
					,
					Message[gate::warning, "Atoms "<>StringRiffle[lostq, ","]<>" are gone to the environment. Return complete depolarising instead on qubits "<>StringRiffle[{q},", "]];
					Subscript[Depol, #][3/4] & /@ {q}
					,
					circ
				]
			];

		(* Vacuum life time limits qubits coherence *)
		t1 = vaclifetime / qubitsnum;
		
		(* standard passive noise *) 
		stdpn[q_, dur_] := 
			{Subscript[Depol, q][0.75 (1 - E ^ (-dur / t1))], Subscript[Deph, q][0.5 (1 - E ^ (-dur / t2))]};	
			
		(* move noise, similar *)				
		movenoise[q_, dur_] := 
			{Subscript[Depol, q][0.75 (1 - E ^ (-dur / t1))], Subscript[Deph, q][0.5 (1 - E ^ (- heatfactor dur / t2))]};
										
		(* coordinate distance measure *)
		distloc[q1_, q2_] := 
			Norm[atomlocs[q1] - atomlocs[q2], 2]  unitlattice;
		
		(* 
		legitimate multi-qubit gates by blockade condition 
		*) 
		blockadecheck[q_List] := If[IntersectingQ[q, Keys @ lossatomlocs], 
			Message[LossAtoms::error, "Some atoms are lost before applying a multi-qubit gate", Return @ False]
			,
			And @@ ((distloc @@ # <= 2 blockaderad)& /@ Subsets[Flatten[q], {2}])
		];
				
		(* 
		legitimate shift move check if the new spots are unoccupied 
		*)
		legshift[q_List, v_] := Module[
			{initlocs = atomlocs},
			initlocs[#] += v & /@ Complement[q, Keys @ lossatomlocs];
			Length @ initlocs === Length @ DeleteDuplicates @ Values @ initlocs
		];
	
	<|
		OptionsUsed -> {
			QubitNum -> qubitsnum, 
			AtomLocations -> atomlocations,
			T2 -> t2,
			VacLifeTime -> vaclifetime,
			RabiFreq -> rabifreq,
			ProbBFRot -> <|10 -> probbfrot10, 01 -> probbfrot01|>,
			UnitLattice -> unitlattice,
			BlockadeRadius -> blockaderad,
			ProbLeakInit -> probleakinit,
			DurInit -> durinit,
			FidMeas -> fidmeas,
			DurMeas -> durmeas,
			DurMove -> durmove,
			ProbLossMeas -> problossmeas,
			ProbLeakCZ -> probleakcz
		}
		,
		DeviceDescription -> ToString[qubitsnum]<>" Rydberg atoms in a "<>ToString@Length[First @ Values @ atomlocs]<>"D lattice."
		,
		NumAccessibleQubits -> qubitsnum
		,
		NumTotalQubits -> qubitsnum
		,
		(** custom keys to access internal variables **)
		LossAtomsProbability :> lossatomprob
		,
		LossAtoms :> lossatomlocs
		,
		AtomLocations :> atomlocs
		,
		BlockadeRadius -> blockaderad
		,
		UnitLattice -> unitlattice
		,	
		(*
		(* re-initialized when invoking InsertCircuitNoise *)
		InitVariables -> Function[
			atomlocs = atomlocations;
			lossatom = <| Table[k -> False, {k, Keys @ atomlocations}] |>;
			lossatomprob = <| Table[k -> 0, {k, Keys @ atomlocations}] |>;
			]
		,
		*)
		(* Aliases are useful for declaring custom operators. At this stage the specification is noise-free (but see later) *)
		Aliases -> {
			Subscript[Init, q_Integer] :> Nothing
			,
			Subscript[SRot, q_Integer][phi_, delta_, tg_] :> Circuit[Subscript[U, q][unitaryDetuning[phi, delta, tg, rabifreq]]]
			,
			Subscript[CZ, p_Integer, q_Integer][phi_] :> Circuit[Subscript[U, p, q][{{1, 0, 0, 0}, {0, E ^ (I phi), 0, 0}, {0, 0, E ^ (I phi), 0}, {0, 0, 0, E ^ (I(2 phi - Pi))}}]]
			,
			Subscript[SWAPLoc, q1_Integer, q2_Integer] :> Nothing
			,
			Subscript[Wait, q__][t_] :> Nothing
			,
			Subscript[ShiftLoc, q__][v_] :> Nothing
			,
			(* multi-qubit gates *)
			Subscript[C, c_Integer][Subscript[Z, t__Integer]] :> Table[Subscript[C, c][Subscript[Z, targ]], {targ, {t}}]
			,
			Subscript[C, c__Integer][Subscript[Z, t_Integer]] :> {Subscript[C, c][Subscript[Z, t]]}
		},
		
		(* Global time: not yet used here *)
		TimeSymbol ->  globaltime
		,
		(* gates rules *)
		Gates -> {
				Subscript[Init, q_Integer] :> 
				<|
				 (* Put the electron back to the atom and reset leak probability *)
					UpdateVariables -> Function[
										lossatomprob[q] = 0;
										KeyDropFrom[lossatomlocs, q];										
									],
					(* perfect init + leakage *)			
					NoisyForm -> {Subscript[Damp, q][1], Subscript[KrausNonTP, q][{{{Sqrt[1 - probleakinit], 0}, {0, 1}}}]}, 
					GateDuration -> durinit
				|>
				,
				Subscript[Wait, q__][dur_] /; (Complement[Flatten @ {q}, qubits] === {}) :>
				 <|
					NoisyForm -> Table[stdpn[j, dur], {j, Flatten @ {q}}],
					GateDuration -> dur
				|>
				,
				Subscript[M, q_Integer] :> 
				<|	
					NoisyForm -> If[MemberQ[Keys @ lossatomlocs, q], 
						Message[M::error, "Atom "<>ToString[q]<>" is gone to the environment, can't measure"],
						{Subscript[Depol, q][Min[1 - fidmeas, 3/4]], Subscript[M, q]}
					], 
					UpdateVariables -> Function[
						lossatomprob[q] = 1 - (1 - lossatomprob[q])(1 - problossmeas);
						(* if the atom is lost, remove it from the location lists *)
						If[\[Not]MemberQ[Keys @ lossatomlocs, q ],
							(* throw a dice to throw the a *)
							If[RandomVariate[BinomialDistribution[1, lossatomprob[q]]]/.{0 -> False, 1 -> True},
								lossatomlocs[q] = atomlocs[q];
								KeyDropFrom[atomlocs, q];
							]
						]
					], 
					GateDuration -> durmeas
				|>
				(* Single rotation by detuning, the native gate. The rest are derivative from SRot *)
				,
				Subscript[SRot, q_Integer][\[Phi]_, \[CapitalDelta]_, tg_] :> 
				<|
					NoisyForm -> circorloss[SRot, {Subscript[SRot, q][\[Phi], \[CapitalDelta], tg], Subscript[Deph, q][0.5( 1 - E^(- tg 0.5 / t2))]}, q],
					GateDuration -> tg
				|>
				,
				Subscript[H, q_Integer] :> 
				<|
					NoisyForm -> circorloss[H, {Subscript[H, q], Subscript[Deph, q][0.5 (1 - E^(-0.5 /(rabifreq t2)))]}, q],
					GateDuration-> \[Pi]/rabifreq
				|>
				,
				Subscript[Rx, q_Integer][\[Theta]_] :> 
				<|
					NoisyForm -> circorloss[Rx, {Subscript[Rx, q][\[Theta]], Subscript[asymBitFlip, q][probbfrot01, probbfrot10], Subscript[Deph, q][0.5(1 - E^(-0.5 Abs[\[Theta]/\[Pi]] / (rabifreq t2)))]}, q],
					GateDuration-> Abs[\[Theta]]/rabifreq
				|>
				,
				Subscript[Ry, q_Integer][\[Theta]_] :> 
				<|
					NoisyForm -> circorloss[Ry, {Subscript[Ry, q][\[Theta]], Subscript[asymBitFlip, q][probbfrot01, probbfrot10], Subscript[Deph, q][0.5(1 - E^(-0.5 Abs[\[Theta]/\[Pi]]/(rabifreq t2)))]}, q],
					GateDuration-> Abs[\[Theta]]/rabifreq
				|>
				,
				Subscript[Rz, q_Integer][\[Theta]_] :>
				<|
					NoisyForm -> circorloss[Rz, {Subscript[Rz, q][\[Theta]], Subscript[asymBitFlip, q][probbfrot01, probbfrot10], Subscript[Deph, q][0.5 (1 - E^(-0.5 Abs[\[Theta]/\[Pi]]/(rabifreq t2)))]}, q],
					GateDuration -> Abs[\[Theta]]/rabifreq
				|>
				, 
				(* moves *)
				Subscript[SWAPLoc, q1_Integer, q2_Integer] :> 
				<|
					UpdateVariables -> Function[
						(* check cases: present atoms, partial loss, complete loss *)
						Which[
							\[Not]MemberQ[Keys @ lossatomlocs, q1] && \[Not]MemberQ[Keys @ lossatomlocs, q2],
							{atomlocs[q1], atomlocs[q2]} = {atomlocs[q2], atomlocs[q1]};
							,
							\[Not]MemberQ[Keys @ lossatomlocs, q1] && MemberQ[Keys @ lossatomlocs, q2],
							atomlocs[q1] = atomlocs[q2];
							Message[SWAPLoc::warning, "Atom "<>ToString[q2]<>" has lost to the environment, only location of atom "<>ToString[q1]<>" changed."]
							,
							MemberQ[Keys @ lossatomlocs, q1] && \[Not]MemberQ[Keys @ lossatomlocs, q2],
							atomlocs[q2] = atomlocs[q1];
							Message[SWAPLoc::warning, "Atom "<>ToString[q1]<>" has lost to the environment, only location of atom "<>ToString[q2]<>" changed."]
							,
							True
							,
							Message[SWAPLoc::warning, "Atom "<>ToString[q1]<>" and "<>ToString[q2]<>" has lost to the environment, no moving on them."]
						]
						
						],
					NoisyForm -> Flatten[circorloss[SWAPLoc, movenoise[#, durmove], #]& /@ {q1, q2}],
					GateDuration -> durmove
				|>
				,
				Subscript[ShiftLoc, q__][v_] /; legshift[{q}, v] :> 
				<|
					UpdateVariables -> Function[																												
						With[{qgone = Intersection[Keys @ lossatomlocs, {q}]}, 
							If[Length @ qgone > 0, 
								Message[ShiftLoc::warning, "Atoms "<>StringRiffle[qgone,", "]<>" are lost to the environment. Do nothing on them."]
							]
						];
						If[KeyExistsQ[atomlocs , #], atomlocs[#] += v] & /@ {q};
					],
					NoisyForm -> Flatten[circorloss[ShiftLoc, movenoise[#, durmove], #]& /@ {q}],
					GateDuration -> durmove
				|>
				,
				(* 
					Multi-qubit gates need all atoms to be present which is checked in the blockadecheck function
				 *)
				Subscript[CZ, p_Integer, q_Integer][\[Phi]_] /; blockadecheck[{p, q}] :> 
				<|
					NoisyForm -> {Subscript[CZ, p, q][\[Phi]], Subscript[KrausNonTP, p, q][{{{1, 0, 0, 0}, {0, Sqrt[1 - probleakcz[01]], 0, 0}, {0, 0, Sqrt[1-probleakcz[01]], 0}, {0, 0, 0, Sqrt[1-probleakcz[11]]}}}]},
					GateDuration -> 4\[Pi]/rabifreq
				|>
				,
				(* The parameterised multi-gates are suspended *)
				Subscript[C, c_][Subscript[Z, t__]] /; blockadecheck[{c, t}] :>
				<|
					NoisyForm -> Join[Table[Subscript[C, c][Subscript[Z, targ]], {targ, {t}}], Subscript[KrausNonTP, #][{{{1, 0}, {0, Sqrt[1 - probleakcz[11]]}}}]& /@ Flatten @ {c, t}],
					GateDuration ->  4\[Pi]/rabifreq
				|>
				,
				Subscript[C, c__][Subscript[Z, t_]] /; blockadecheck[{c, t}] :> 
				<|
					NoisyForm -> Join[{Subscript[C, c][Subscript[Z, t]]}, Subscript[KrausNonTP, #][{{{1, 0}, {0, Sqrt[1 - probleakcz[11]]}}}]& /@ Flatten @ {c, t}],
					GateDuration -> 4\[Pi]/rabifreq
				|>
			}
		,
		DurationSymbol -> \[CapitalDelta]t,
			Qubits -> { 
				q_ :> <| 
					PassiveNoise -> stdpn[q, \[CapitalDelta]t]
				|>
			}
		|>
		]
	]

		
	
	(* DEVICE_SILICONDELFT *)
	
	SiliconDelft[OptionsPattern[]] := With[
	{
		qubitnum = OptionValue @ QubitNum,
		efsinglexy = OptionValue @ EFSingleXY,
		efcz = OptionValue @ EFCZ,
		t1 = OptionValue @ T1,
		t2 = OptionValue @ T2,
		rabifreq = OptionValue @ RabiFreq,
		qubitfreq = OptionValue @ QubitFreq,
		freqcz = OptionValue @ FreqCZ,
		fidcz = OptionValue @ FidCZ,				
		fidsinglexy = OptionValue @ FidSingleXY,
		fidread = OptionValue @ FidRead,
		fidcrot = OptionValue @ FidCRot,
		durread = OptionValue @ DurRead,
		freqcrot = OptionValue @ FreqCRot,
		(* True/False*)
		exchangerotoff = OptionValue @ ExchangeRotOff,
		offresonantrabi = OptionValue @ OffResonantRabi,
		stdpassivenoise = OptionValue @ StdPassiveNoise,
		(*a matrix or a boolean*)
		exchangeroton = OptionValue @ ExchangeRotOn,
		efsinglexy = OptionValue @ EFSingleXY,
		efcz = OptionValue @ EFCZ	
	}
	,
	
	(* 
		assertions and formatter (TODO)  
	*)
	
	(* even number of qubits *)
	If[\[Not]EvenQ[qubitnum], 
		Message[QubitNum::"The number of qubits must be even."; Return @ $Failed]];
		
	
	Module[
	{er1xy, ercz, qubits, deltaT, passivenoisecirc, offresrabi, stdpn, exczon, sroterr, g2, ndeph, ndepol, a, sq, eq, first, second, ccrot, edgeq},
		
		(* conveniently defined stuff *)
		qubits = Range[0, qubitnum - 1];
		
		(* error parameters *)
		er1xy = fid2DepolDeph[#, efsinglexy, 1, FidSingleXY, True]& /@ fidsinglexy;
		ercz = fid2DepolDeph[#, efcz, 2, FidCZ, True]& /@ fidcz;

		a = qubitnum; (* index of the ancilla qubit for parity measurement purpose  *)
		sq = 0; (* index of start edge*)
		eq = qubitnum - 1; (* index of ending edge *)
		g2 = False; (* indicate if a 2-qubit gate is actively applied *)
		
		(*first half, second half of the qubits*)
		{first, second} = Partition[qubits, qubitnum / 2];
		
		(* Legal Subscript[CROT, c,t] operations *)
		ccrot = Join[Partition[Reverse @ first, 2, 1], Partition[second, 2, 1]];
		
		(* edge qubits *)
		edgeq = {first[[;;2]], Reverse @ first[[;;2]], second[[-2;;]], Reverse @ second[[-2;;]]};
		
		(* Normalise the numbers to be within the correct range of error parameters*)		
		ndeph[num_] := Min[num, 0.5];
		ndepol[num_] := Min[num, 0.75];	
		
		(* standard passive noise *)
		stdpn[q_, dur_] := If[stdpassivenoise,
			{Subscript[Deph, q][.5(1 - Exp[-N[dur/t2[q], $MachinePrecision]])], Subscript[Depol, q][0.75(1 - Exp[-dur/t1[q]])]}
			,
			{}];
			
		(* the entire passive noise *)
		passivenoisecirc[q_Integer, g2_, dur_] := Flatten @ {
		If[\[Not]g2 && q < qubitnum - 1 && AssociationQ @ exchangerotoff, (**)
			Subscript[C, q][Subscript[Rz, q+1][(dur/\[Pi])*exchangerotoff[q]]]
			,
			{}
		], stdpn[q, dur]};
		
		(* cross-talk single rotation noise due to detuning *)
		offresrabi[q_, \[Theta]_]:= If[offresonantrabi, 
			Table[Subscript[U, j][offresonantRabi[rabifreq[q], qubitfreq[j] - qubitfreq[q], Abs[\[Theta]]]],{j, Delete[qubits, q+1]}] 
			,
			{}];
		
		(*Exchange rotation C-Rz[j] interaction when CZ gate on*)
		exczon[targ_] := If[ListQ @ exchangeroton,
			Subscript[C, #-1][Subscript[Rz, #][exchangeroton[[targ, #]]]]& /@ Delete[Range[qubitnum-1],targ]
			,
			{}];
		
		(*Errors on single rotations*)
		sroterr[q_,\[Theta]_] := Flatten @ {
			Subscript[Depol, q][ndepol[er1xy[q][[1]] Abs[\[Theta]/\[Pi]]]], Subscript[Deph, q][ndeph[er1xy[q][[2]] Abs[\[Theta]/\[Pi]]]], offresrabi[q, \[Theta]]};						
		
		<|
		(* one hidden qubit for measurement ancilla *)
		DeviceType -> "SiliconDelft"
		,
		OptionsUsed -> {
			QubitNum -> qubitnum,
			
		efsinglexy = OptionValue @ EFSingleXY,
		efcz = OptionValue @ EFCZ,
		t1 = OptionValue @ T1,
		t2 = OptionValue @ T2,
		rabifreq = OptionValue @ RabiFreq,
		qubitfreq = OptionValue @ QubitFreq,
		freqcz = OptionValue @ FreqCZ,
		fidcz = OptionValue @ FidCZ,				
		fidsinglexy = OptionValue @ FidSingleXY,
		fidread = OptionValue @ FidRead,
		fidcrot = OptionValue @ FidCRot,
		durread = OptionValue @ DurRead,
		freqcrot = OptionValue @ FreqCRot,
		(* True/False*)
		exchangerotoff = OptionValue @ ExchangeRotOff,
		offresonantrabi = OptionValue @ OffResonantRabi,
		stdpassivenoise = OptionValue @ StdPassiveNoise,
		(*a matrix or a boolean*)
		exchangeroton = OptionValue @ ExchangeRotOn,
		efsinglexy = OptionValue @ EFSingleXY,
		efcz = OptionValue @ EFCZ		
		}
		,
		DeviceDescription -> "Silicon spin device Delft-inspited with "<>ToString[qubitnum]<>"-qubits arranged as a linear array with nearest-neighbor connectivity. One extra qubit is used as an ancilla to simulate measurement."
		,
		NumAccessibleQubits -> qubitnum
		,
		NumTotalQubits -> qubitnum + 1
		,		
		Aliases -> {
			Subscript[Wait, q__][t_] :> Sequence@@{},
			Subscript[CRot, q0_,q1_]:> Subscript[C, q0][Subscript[X, q1]],
			Subscript[MPar, q0_,q1_ ]:> Sequence@@{Subscript[Damp, a][1],Subscript[X, a],Subscript[H, a],Subscript[C, a][Subscript[Z, q0]],Subscript[C, a][Subscript[Z, q1]],Subscript[H, a],Subscript[M, a],
			Subscript[Kraus, q0,q1][{{{1,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,1,0},{0,0,0,0},{0,0,0,0}},
			{{0,0,0,0},{0,1,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}}}]}
			}
			,	
		Gates ->{
			Subscript[Wait, q__][t_]:><|
				NoisyForm-> Flatten@Table[passivenoisecirc[i,False,t],{i,Flatten@{q}}],
				GateDuration->t,
				UpdateVariables->Function[g2=False]
			|>,
		(* Parity measurement only *)
			Subscript[MPar, q0_,q1_]/; MemberQ[edgeq,{q0,q1}] :><|
				NoisyForm-> {Subscript[Kraus, q0,q1][bitFlip2[fidread]],Subscript[MeasP, q0,q1],Subscript[Kraus, q0,q1][bitFlip2[fidread]]},
				GateDuration->durread,
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
			Subscript[C, p_][Subscript[Z, q_]]/; Abs[q-p]==1  :><|
				(*The last bit undo the exchange in the passive noise *)
				NoisyForm->{Subscript[C, p][Subscript[Z, q]],Subscript[Depol, p,q][ercz[p][[1]]],Subscript[Deph, p,q][ercz[p][[2]]],Sequence@@exczon[q]}, 
				GateDuration->\[Pi]/freqcz[p],
				UpdateVariables->Function[g2=True]
			|>,
			
			Subscript[C, p_][Subscript[Ph, q_][\[Theta]_]]/; Abs[q-p]==1  :><|
				(*The last bit undo the exchange in the passive noise *)
				NoisyForm->{Subscript[C, p][Subscript[Ph, q][\[Theta]]],Subscript[Depol, p,q][Min[ercz[p][[1]]Abs[\[Theta]/\[Pi]],15/16]],Subscript[Deph, p,q][Min[ercz[p][[2]]*Abs[\[Theta]/\[Pi]],3/4]],Sequence@@exczon[q]}, 
				GateDuration->Abs[\[Theta]]/freqcz[p],
				UpdateVariables->Function[g2=True]
			|>,
			Subscript[CRot, c_,t_]/; MemberQ[ccrot,{c,t}]  :><|
				NoisyForm->{Subscript[CRot, c,t],Subscript[Depol, c,t][Min[1-fidcrot,15/16]]}, 
				GateDuration->1/freqcrot,
				UpdateVariables->Function[g2=False]
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
	erent=entFid2DepolDeph[OptionValue[FidEnt],OptionValue[EFEnt],FidEnt]
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
	stdp[node_,t_,q_]:={Subscript[Damp, qmap[node][q]][(1-E^(-t/t1[node]))],Subscript[Deph, qmap[node][q]][.5(1-Exp[-(t/t2s[node])^2])]};
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

(** CIRCUIT_CONSTRUCTION **)
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

GenerateOptionTable[options_,columnwidth_:{"3cm","3cm","10cm"}]:=Module[{vars, contents, header,content,contenf,contenl,fvalues,finfo,cw1,cw2,cw3},
	vars=Keys@options;
	{cw1,cw2,cw3}=columnwidth;
	header=StringForm["\\begin{tabular}{p{``}p{``}p{``}}\n\\toprule",cw1,cw2,cw3];
	header=StringRiffle[{header,"\\textbf{Variable}&\\textbf{Value}&\\textbf{Description}\\\\ \n\\midrule\n"}];
	fvalues=List@@#&/@Values@options;
	finfo=ToString@Information[#]&/@Keys[options];
	finfo=Quiet[StringReplace[#,
	{"Rx"->"$Rx$", "Ry"->"$Ry$", "Rz"->"$Rz$","T1"->"$T_1$","T2s"->"$T_2^s$","T2"->"$T_2$"}]&/@finfo;
	content=Transpose@{Keys@options,fvalues,finfo}];
	contenl=StringRiffle[#,"&"]&/@content;
	contents=StringRiffle[contenl,"\\\\ \n"];
	header<>contents<>"\\\\ \n \\bottomrule\n\\end{tabular}"
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

Needs["VQD`CustomGates`"]
