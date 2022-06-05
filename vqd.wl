(* ::Package:: *)

BeginPackage["VQD`"];


Needs["QuEST`"];


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


Options[SiliconDelft]={

};


BeginPackage["`ParameterDevices`"];


(* Global variables *)
DrawIons::usage="Draw the current string of ions";
Init::usage="Initialise qubit to state |0>";
MSCrossTalk::usage="(entangling OR starkshift) The crosstalk model in applying M\[OSlash]lmer\[Dash]S\[OSlash]rensen gate";
NIons::usage="The total number of ions";
DurSingleGate::usage="Duration of single rotation of \[Pi]";
DurTwoGate::usage="Duration of two qubit gates with rotation of \[Pi]";
DurMeas::usage="Duration of measurement";
DurInit::usage="Duration of initialisation";
DurShuffle::usage="Duration to shuffle location of ions";
entangling::usage="Crosstalk error model pre equation (4) on applying the XX MS-gate";
ErrCT::usage="Error coefficient of the crosstalk with entanglement model";
ErrSS::usage="Error coefficient of the crosstalk with stark shift model";
FidSingleXY::usage="Fidelity of single Rx[\[Theta]] and Ry[\[Theta]] rotations obtained by random benchmarking";
FidMeas::usage="Fidelity of measurement";
FidInit::usage="Fidelity of qubit initialisation";
starkshift::usage="Crosstalk error model using stark shift on  applying the XX MS-gate";
T1::usage="T1 duration in \[Mu]s. Exponential decay time for the state to be complete mixed.";
T2::usage="T2 duration in \[Mu]s. Exponential decay time for the state to be classical.";


(* Custom gates *)
SWAPLoc::usage="Swap the spatial locations of two qubits";
Wait::usage="Wait gate, doing nothing";


EndPackage[];


(*All definitions of modules*)


Begin["`Private`"];


End[];


EndPackage[];


Needs["VQD`ParameterDevices`"]
