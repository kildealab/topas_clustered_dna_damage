# clustered_dna_damage_indirect 

clustered_dna_damage_indirect is an extension of topas_clustered_dna_damage to include indirect action modeling and damage scoring.

Updates to topas_clustered_dna_damage:

- ScoreClusteredDNADamage.cc (DNA damage scorer):
  - Added: simulation of indirect action events.
  - Added: indirect damage scoring.
  - Added: simulation constraints in the chemical stage:
    - ·OH radical tracks were terminated after an indirect action event (whether or not DNA damage was inflicted.
    - Radical tracks (·OH, e− aq, and H· specifically) were terminated immediately upon diffusion into a histone volume.
    - The generation of reactive chemical species was not permitted inside of DNA and histone volumes because these volumes are not made of liquid water in reality.
  - Added: user-modifiable simulation parameters:
    - Toggle to score direct damage.
    - Toggle to score indirect damage.
    - Damage probabilities for molecule-base and molecule-backbone interactions.
    - Toggle for histone scavenging.
    - Molecule species scavenged by the histone volumes.
    - Molecule species scavenged by the DNA volumes.
  - Updated: multithreading feature to support indirect action and damage scoring.
  - Updated: DNA damage clustering algorithm to account for indirect and hybrid lesions.
  
- VoxelizedNuclearDNA.cc (DNA model):
  - Added: unique identification of histone volumes via their composing material.
