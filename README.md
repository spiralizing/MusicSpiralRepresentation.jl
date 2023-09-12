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

### To load MIDI (CSV) files

You can also load MIDI files indirectly by converting them to CSV files first with [midicsv](https://www.fourmilab.ch/webtools/midicsv/) in a terminal:

```
$ midicsv name_of_file.mid > name_of_file.csv 
```

and to load it in Julia

```julia
using DelimitedFiles

csv_piece = readdlm("name_of_file.csv",',')

```

### Creating a DataFrame from MusicXML and CSV files

You can create a `DataFrame` containing the information of the music score/piece with `get_xml_df` or `get_csv_df`, where the resulting data frame returns the format:

| `measure`| `Time signature`|`start` | `end` | `duration` | `pitch`| 

for xml files:

```julia
import MusicSpiralRepresentation
const msr = MusicSpiralRepresentation
using PyCall
m21 = pyimport("music21")

mxl_piece = m21.converter.parse("name_of_file.mxl")
piece_df = msr.get_xml_df(mxl_piece)
```

for csv files:

```julia
import MusicSpiralRepresentation
const msr = MusicSpiralRepresentation
using DelimitedFiles

csv_piece = readdlm("name_of_file.csv",',')
piece_df = msr.get_csv_df(csv_piece)
```
