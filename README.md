# TOPAS_Clustered_DNA_Damage

![Logo](https://github.com/McGillMedPhys/clustered_dna_damage/blob/dev/repository_logo_figure.svg)

This repository contains a TOPAS-nBio application that can be used to simulate clustered DNA damage due to direct and indirect action of ionizing radiation.

v1: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5090104.svg)](https://doi.org/10.5281/zenodo.5090104)

v2:

## Table of Contents

* [Authors](#authors)
* [Features](#features)
* [Description](#description)
* [Dependencies](#dependencies)
* [Installation](#installation)
* [Instructions](#instructions)
* [Output](#output)
* [License](#license)
* [Component Details](#component-details)
* [Changes from Last Version](#changes-from-last-version)

## Authors

Logan Montgomery, Christopher M Lund, James Manalad, Anthony Landry, John Kildea

Contact email: logan.montgomery@mail.mcgill.ca, james.manalad@mail.mcgill.ca

## Features

* Complete TOPAS parameter file required to run simulations.
* Full human nuclear DNA model (implemented as a custom geometry component).
* Algorithm to record clustered DNA damage (implemented as a custom scorer).
* Physics constructor (implemented as a custom physics module).
* Energy spectra and relative dose data files for secondary particles produced by neutrons and x-rays in human tissue.
* All code is thoroughly documented.

## Description

* This application is intended to be used to simulate the induction of clustered DNA damage in a human nucleus.
* We developed this application to compare neutron-induced direct and indirect clustered DNA damage with x-ray induced DNA damage in order to invesigate the energy dependence of neutron RBE.
* Specifically, the application produces yields of the following DNA damage:
    1. Single strand breaks (SSBs)
    2. Base lesions
    3. Double strand breaks (DSBs)
    4. Complex DSB clusters (clusters containing at least 1 DSB).
    5. Non-DSB clusters (clusters that don't contain any DSBs).
* Most simulation parameters can be modified using the included [parameter file](https://github.com/McGillMedPhys/topas_clustered_dna_damage/blob/indirect/DNAParameters.txt).
* Details about each component of this application are provided [below](#component-details).

## Dependencies

* TOPAS v3.6.1
* TOPAS-nBio 1.0

**Note**: This application was developed on Ubuntu 20.04.2.

## Installation

1. Download the latest version from the [releases page](https://github.com/McGillMedPhys/clustered_dna_damage/releases).
2. Install the [dependencies](#dependencies).
3. Install TOPAS_Clustered_DNA_Damage as any other TOPAS extension as per the [instructions provided by TOPAS](https://sites.google.com/a/topasmc.org/home/home).
    1. Place this repository in your `topas_extensions` directory.
    2. Recompile TOPAS, e.g:
        * `cd /path/to/topas`
        * `cmake -DTOPAS_EXTENSIONS_DIR=/path/to/topas_extensions`
        * `make`

## Instructions

1. Enter desired settings for the application by editing the parameter file (`DNAParameters.txt`)
2. Run the application (`topas DNAParameters.txt`)

## Output

| File | Description |
| ----------- | ----------- |
| damage_yields.phsp | Yields of [five types of DNA damage](#description) stratified according to their damage cause: direct action, indirect action, or both (hybrid)|
| run_summary.csv | Details about the simulation run |
| data_comp_dsb.csv | Cluster properties of every recorded complex DSB cluster |
| data_non_dsb.csv | Cluster properties of every recorded non-DSB cluster  |

## License

* This project is provided under the MIT license. See the [LICENSE file](LICENSE) for more info.
* When using any component of this application, please be sure to cite our papers:
    * Montgomery L, Lund CM, Landry A, Kildea J (2021). Towards the characterization of neutron carcinogenesis through direct action simulations of clustered DNA damage. <em>Phys Med Biol</em> 66(20); 205011.
        * DOI: [https://doi.org/10.1088/1361-6560/ac2998](https://doi.org/10.1088/1361-6560/ac2998)
    * Manalad J, Montgomery L, Kildea J (2022). (coming soon)
        * DOI: (coming soon)

## Component details

### Nuclear DNA model
* Source code file is located [here](https://github.com/McGillMedPhys/clustered_dna_damage/blob/master/geometry/VoxelizedNuclearDNA.cc).
* Full human nuclear DNA model containing ~6.3 Gbp.
* Cubic shape constructed using voxels.
* Each voxel contains 20 chromatin fibres.
* Every fibre contains 18,000 DNA base pairs.
* Nucleus is enclosed in a spherical cell volume (fibroblast model).

### Clustered DNA damage scorer
* Source code file is located [here](https://github.com/McGillMedPhys/clustered_dna_damage/blob/master/scoring/ScoreClusteredDNADamage.cc).
* Simulates direct and indirect prompt DNA damage.
* During the chemical stage:
    * All radical tracks generated inside DNA and histone volumes are immediately terminated.
    * DNA and histone volumes can "scavenge" (terminate) radiolytic species.
* Records the five types of DNA damage [mentioned above](#description) and their respective damage-inducing action.
* Damage definitions (separation distances, energy thresholds, indirect damage probabilities) can be modified in the parameter file as shown [here](https://github.com/McGillMedPhys/topas_clustered_dna_damage/blob/indirect/supportFiles/DNADamageParameters.txt).
* Other user-modifiable simulation parameters:
    * Toggles to score direct and indirect damage, and histone scavenging.
    * Molecule species scavenged by the DNA and histone volumes.
* Default behaviour is to terminate simulation after a fixed number of histories.
    * Can alternatively terminate simulation after a certain dose deposition in the nucleus.
* Supports multithreading.
* Default parameter values related to indirect action and the chemical stage are described [below](#changes-from-last-version).

### Physics module
* Source code file is located [here](https://github.com/McGillMedPhys/clustered_dna_damage/blob/master/physics/G4EmDNAPhysics_option2and4.cc).
* Combines the GEANT4-DNA physics constructors: `G4EmDNAPhysics_option2` and `G4EmDNAPhysics_option4`.
* Physics models from `G4EmDNAPhysics_option4` for electrons between 10 eV and 10 keV.
* Physics models from `G4EmDNAPhysics_option2` for electrons between 10 keV and 1 MeV.

### Secondary particle data files
* In a previous study, we evaluated the energy spectra and relative dose contributions of secondary particles produced by neutrons & 250 keV x-rays in human tissue.
* For details, see our paper:
    * Lund CM, Famulari G, Montgomery L, Kildea J (2020). A microdosimetric analysis of the interactions of mono-energetic neutrons with human tissue. <em>Physica Medica</em> 73; 29-42.
        * DOI: [https://doi.org/10.1016/j.ejmp.2020.04.001](https://doi.org/10.1016/j.ejmp.2020.04.001)
* These data are included as TOPAS parameter files in this repository.
    * Spectra are located [here](https://github.com/McGillMedPhys/clustered_dna_damage/tree/master/spectra).
    * Relative dose values are located [here](https://github.com/McGillMedPhys/clustered_dna_damage/tree/master/relative_doses).
* Naming convention of these files:
    * e.g. `spectrum_n1MeV_inner_proton.txt`
        * `n1MeV`: initial 1 MeV neutrons.
        * `inner`: irradiated the innermost scoring volume in human tissue.
        * `proton`: protons produced as secondary particles.
 * These files can be referenced in the main parameter file `DNAParameters.txt` to irradiate the nuclear DNA model.

## Changes from last version

### Nuclear DNA model:
* Unique identification of histone volumes via their composing material was added.

### Clustered DNA damage scorer:
* Simulation of indirect action events and indirect damage scoring using the model described in:
    * Zhu H _et al_. (2020). Cellular response to proton irradiation: a simulation study with TOPAS-nBio. <em>Radiation Research</em> 194; 9-21.
        * DOI: [https://doi.org/10.1667/rr15531.1](https://doi.org/10.1667/rr15531.1)
* Constraints simulated by default during the chemical stage:
    * All radical tracks generated inside DNA and histone volumes are immediately terminated.
    * ·OH radical tracks are terminated after an indirect action event (whether or not DNA damage was inflicted).
    * Radical tracks (·OH, e<sup>-</sup><sub>aq</sub>, and H· specifically) are terminated immediately upon diffusion into a histone volume.
* By default, only ·OH radicals can damage DNA volumes with a damage probability of 40%.
    * The damage probabilities of other radiolytic species with backbone or nitrogenous base volumes can be modified via the parameter file.
* Other user-modifiable simulation parameters:
    * Toggle to score direct damage.
    * Toggle to score indirect damage.
    * Toggle for histone scavenging.
    * Molecule species scavenged by the DNA volumes.
    * Molecule species scavenged by the histone volumes.
* The DNA damage clustering algorithm was updated to account for indirect and hybrid lesions.
* Multithreading support for indirect action simulations to decrease simulation time.
