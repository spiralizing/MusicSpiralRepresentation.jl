# MusicSpiralRepresentation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://spiralizing.github.io/MusicSpiralRepresentation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://spiralizing.github.io/MusicSpiralRepresentation.jl/dev/)
[![Build Status](https://github.com/spiralizing/MusicSpiralRepresentation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/spiralizing/MusicSpiralRepresentation.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/spiralizing/MusicSpiralRepresentation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/spiralizing/MusicSpiralRepresentation.jl)

Julia Package with functions to compute the Center of Effect with the Spiral Representation as model for tonality developed by Elaine Chew in her book *Mathematical and Computational Modeling of Tonality*.

## Install

```julia-repl
pkg> add https://github.com/spiralizing/MusicSpiralRepresentation.jl

```

## Usage
For a bit of explanation of what this package is about and how to use some of its functions please see [THIS EXAMPLE](https://spiralizing.github.io/DSEntries/CenterOfEffect/)

### To load MusicXML files 
**You need to have PyCall installed in your Julia environment**

```julia
import MusicSpiralRepresentation
const msr = MusicSpiralRepresentation

using PyCall
m21 = pyimport("music21")

mxl_piece = m21.converter.parse("name_of_file.mxl")
```
