module MusicSpiralRepresentation
using DataFrames
using Distances
using EzXML
using LinearAlgebra
using MIDI
#using PyCall #uncomment this line if you want 
# to use the python library music21 as fallback.
using Statistics
using ZipFile
#
export pitch_names, nminor_keys, nMajor_keys, cf_notes
export nall_keys, midi_notes
export get_center_effect, get_distance_to_keys, get_distance_to_chords, get_piece_by_measure, get_csv_df, get_xml_df
export get_piece_by_measure_mid, get_mid_df
export get_piece_by_measure_xml
export get_rank_freq

include("Utils.jl")
include("SpiralDefinitions.jl")
include("CenterOfEffect.jl")

end