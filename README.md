<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/Strong-Lab/Phybrid">
    <img src="images/logo.png" alt="Logo" width="400" height="100">
  </a>
  <p align="center">
    A viral identification tool using machine learning with nucleotide and protein features
    <br />
    <a href="https://github.com/Strong-Lab/Phybrid/issues">Report Bug</a>
    ·
    <a href="https://github.com/Strong-Lab/Phybrid/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-Phybrid">About Phybrid</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About Phybrid

This project was created to identify viral contigs in metagenomics. Add details about project and images

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple example steps.

### Prerequisites

Phybrid requires Python 3 and the following libraries (if installling through pip, libraries are automatically install)
* pandas
* scikit-learn
* biopython


### Installation

Phybrid can only be forked from this repository. Once forked, enter into the files and compile the kmer counting program. 

```
## git download here

cd Phybrid/scripts
tar xvf kmer-counter-master.zip
cd kmer-counter-master
make
```

<!-- USAGE EXAMPLES -->
## Usage
Phybrid works as a python script. Once install via pip, Phybrid the command can be accessed. To get the help screen type:
```
Phybrid -h
```

The paramters of Phybrid are:
* -i: Input Fasta \[required]
* -o: Output Directory \[optional]


### Running Phybrid 

#### Phybrid without a sequencing file and renaming the output
```
cd Phybrid
scripts/Phybrid.py -i data/Test/Viral_contigs.fasta -o data/Test/Output
```


<!-- ROADMAP -->
## Roadmap

<p align="center">
    Current Version: 0.0.1
</p>

Improvements to be made:
- Reduce feature space to allow for smaller file processing
- Grid search hyper parameters for models
- Build into python package


See the [open issues](https://github.com/othneildrew/Best-README-Template/issues) for a list of proposed features (and known issues).


<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- CONTACT -->
## Contact

Cody Glickman - [@glickman_Cody](https://twitter.com/glickman_cody) - glickman.cody@gmail.com

Project Link: [https://github.com/Strong-Lab/Phybrid](https://github.com/Strong-Lab/Phybrid)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* James Costello
* Michael Strong
* Jo Hendrix





<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/Strong-Lab/Phybrid.svg?style=for-the-badge
[contributors-url]: https://github.com/ontributors/Strong-Lab/Phybrid/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/Strong-Lab/Phybrid.svg?style=for-the-badge
[forks-url]: https://github.com/Strong-Lab/Phybrid/network/members
[stars-shield]: https://img.shields.io/github/stars/Strong-Lab/Phybrid.svg?style=for-the-badge
[stars-url]: https://github.com/Strong-Lab/Phybrid/stargazers
[issues-shield]: https://img.shields.io/github/issues/Strong-Lab/Phybrid.svg?style=for-the-badge
[issues-url]: https://github.com/Strong-Lab/Phybrid/issues


