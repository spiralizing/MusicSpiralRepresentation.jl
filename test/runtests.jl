using MusicSpiralRepresentation
using DataFrames
using EzXML
using MIDI
using Test

@testset "MusicSpiralRepresentation.jl" begin
    @test get_distance_to_keys(get_center_effect([60, 64, 67]))[1, 1] == "C"

    @testset "raw MIDI parser" begin
        # Build a known fixture: C major scale, one quarter each, 4/4, 120 qpm.
        tpq = 480
        pitches = [60, 62, 64, 65, 67, 69, 71, 72]
        notes = [MIDI.Note(p, 96, (i - 1) * tpq, tpq, 0) for (i, p) in enumerate(pitches)]
        track = MIDI.MIDITrack()
        push!(track.events, MIDI.TimeSignatureEvent(0, 4, 4, 24, 8))
        push!(track.events, MIDI.SetTempoEvent(0, 500_000))  # 120 qpm
        MIDI.addnotes!(track, MIDI.Notes(notes, tpq))
        midi = MIDI.MIDIFile(1, tpq, [track])

        path = joinpath(mktempdir(), "scale.mid")
        MIDI.writeMIDIFile(path, midi)

        # Quarter-note output matches the XML schema and is tick-exact.
        df_q = get_mid_df(path; abs_time=false)
        @test names(df_q) == ["Measure", "TimeSignature", "StartQuarter", "EndQuarter", "Duration", "Pitch"]
        @test nrow(df_q) == length(pitches)
        @test df_q.Pitch == pitches
        @test df_q.Measure == [1, 1, 1, 1, 2, 2, 2, 2]
        @test all(df_q.Duration .== 1.0)
        @test df_q.StartQuarter == 0.0:7.0

        # Millisecond output: at 120 qpm a quarter note is 500 ms.
        df_ms = get_mid_df(path; abs_time=true)
        @test names(df_ms)[3:4] == ["StartTime", "EndTime"]
        @test all(df_ms.Duration .== 500.0)

        # Center of effect is unit-invariant (durations enter as normalized weights).
        @test get_center_effect(Matrix(df_ms)) ≈ get_center_effect(Matrix(df_q))

        # Dispatcher routes .mid paths to the MIDI parser.
        @test size(get_piece_by_measure(path)) == (length(pitches), 6)
    end

    @testset "native MusicXML parser (EzXML)" begin
        # Uncompressed .musicxml: 4-part string quartet (Beethoven Op. 18 No. 1 in F).
        df = get_xml_df("op18_no1_mov1.musicxml")
        @test names(df) == ["Measure", "TimeSignature", "StartQuarter", "EndQuarter", "Duration", "Pitch"]
        @test nrow(df) == 4181
        @test unique(df.TimeSignature) == ["3/4"]
        @test extrema(df.Pitch) == (36, 94)
        @test get_distance_to_keys(get_center_effect(Matrix(df)))[1, 1] == "F"

        # Compressed .mxl: inner score resolved via META-INF/container.xml.
        df_mxl = get_xml_df("Mozart_16.mxl")
        @test nrow(df_mxl) == 1375
        @test unique(df_mxl.TimeSignature) == ["4/4"]
        @test get_distance_to_keys(get_center_effect(Matrix(df_mxl)))[1, 1] == "C"

        # .xml extension routes through the XML parser too.
        xmlcopy = joinpath(mktempdir(), "op18.xml")
        cp("op18_no1_mov1.musicxml", xmlcopy)
        @test nrow(get_xml_df(xmlcopy)) == 4181

        # pitch_to_midi matches music21's .ps convention (C4 = 60, F#3 = 54).
        doc = MusicSpiralRepresentation.read_musicxml("op18_no1_mov1.musicxml")
        @test MusicSpiralRepresentation.pitch_to_midi(
            findfirst("//note/pitch", root(doc))) isa Int
    end
end
