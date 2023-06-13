#Vitual Quantum Devices (VQD)

**VQD** is a **Mathematica** package contains a collection of virtual quantum devices that are based on or inspired from actual quantum devices. **VQD** is built atop [**QuESTlink**](https://github.com/QTechTheory/QuESTlink), a Mathematica extension of [**QuEST**](https://github.com/QuEST-Kit/QuEST) &mdash; a powerful open-source emulator of quantum computers developed in **C** and **C++**.   [**QuESTlink**](https://github.com/QTechTheory/QuESTlink) combines **Mathematica**’s powerful symbolic operations with [**QuEST**](https://github.com/QuEST-Kit/QuEST)’s high-performance backend, enabling virtual devices to be highly configurable through an intuitive interface, and able to leverage powerful visualisation facilities, without compromising performance.

**VQD provides five families of virtual devices:**
1. Multi-nodes ion traps
2. Nitrogen-vacancy-center  (NV-center) diamond 
3. Rydberg neutral atoms
4. Quantum dots in Si/SiGe 
5. Superconducting transmon qubits.


Details on each quantum device's architecture and error model are available in this **paper**. 


##Getting started

Every quantum device is unique. The virtual device is designed to reflect a close approximation to the physical reality by only providing access to physically possible native operations. Get to know them easily using our straightforward guides in the ``devices`` folder. 

Simply run the notebook of the device of your interest and configure the numbers or characteristics as you wish.
1. **TrappedIonOxford.nb** for the multi-nodes ion traps
2. **NVCenterDelft.nb** for the NV-center qubits
3. **RydbergHub.nb** for the Rydberg neutral atoms
4. **SiliconDelft.nb** for the Silicon spin qubit
5. **SuperconductingHub.nb** for the superconducting qubits.


Each Mathematica notebook helps you understand these devices and their special features. It also gives access to obtain the results shown in our **paper**.
To ensure the correct results when running the notebook, check the files in ``doc`` with the corresponding pdf name.

**Note**: the results shown in **paper** require extensive simulations, which are accessible in folder ``supplement``. Some codes are in text format, to prevent unecessary accidental expensive simulation. To use the code, simply convert the text into code or input.

###General usage

A virtual device contains descriptions  of the devices which comprises the architecture, connectivity, errors, gates frequency, and the device state.  

####Using parallelism, multithreading 

####Navigating through VQD
To use parallelisation
SetEnvironment[] 
sets the system variable, such as:
	the number of threads in multithreading ("OMP_NUM_THREADS"->"#thread") 
	the gpu address ("CUDA_VISIBLE_DEVICES"->"#id") 
This works only for the self-compiled binary 

Example of a VQD setup

UsedOptions

error ratios : sum to 1 or zero (off)