#import PyCall
#m21 = PyCall.pyimport("music21")
"""
    get_piece_by_measure_csv(s; qdiv=32, abs_time=true)

    Returns a vector of two dimensional arrays of the basic information of the piece for analyzing its properties,
    input "s" should be a CSV table from a CSV file converted from a MIDI file with the midicsv script available for free
    in: https://www.fourmilab.ch/webtools/midicsv/ .

    The first output consist of a vector which components represent each measure of the music piece,
    each measure contains four columns; Beat: number of beat where the note starts (beat or halfbeat)
    Measure: the fraction of the measure (n/d), Duration: the duration of the note in midiclocks, if abs_time=true
    the output returns the time in ms where the note started and ended before the duration of the note,
    Pitch: the pitch of the note represented in MIDI notation (0-127).

    The second output corresponds to the total number of notes in the piece and the number of notes that
    fall outside the threshold given (1/qdiv), most of the time if notes fall outside this threshold is because
    the MIDI was sequenced (recorded) instead of generated with a music score software.
"""
function get_piece_by_measure_csv(s; qdiv=32, abs_time=true)
    frac = s[findall(x -> x == " Time_signature", s[:, 3]), 2:5]  #this is the measure of the piece in fractional way a / b
    if length(frac) > 4
        bot = map(x -> 2^x, frac[:, 4])
        sym_measure = map((x, y) -> join([x y], "/"), frac[:, 3], bot)
        measure = map((x, y) -> x / 2^y, frac[:, 3], frac[:, 4])
        meas_time = frac[:, 1]
        numer = frac[:, 3]
    else
        measure = frac[3] / 2^frac[4]
        numer = []
        push!(numer, frac[3])
        sym_measure = join([frac[3] 2^frac[4]], "/")
    end

    q_t = s[1, 6] #the time in ms for a quarter of note 1/4
    th = q_t / qdiv #64nd of tolerance with qdiv=32
    w_s = measure * q_t / 0.25 #this is the miliseconds of a measure
    smdiv = 32 #this is the threshold for the tolerance in the beats
    nv = s[findlast(x -> x == " Note_on_c", s[:, 3]), 1] #- 1 #estimates how many voices are in the midi
    voces = Array{Matrix}(undef, nv) #initialze an array
    subdiv = map((x, y) -> x / y, w_s, numer)
    notes_measure = []
    for i = 1:length(measure)
        push!(notes_measure, cumsum([subdiv[i] for j = 1:numer[i]]))
    end
    n_off = findall(x -> x == " Note_off_c", s[:, 3])
    n_on = findall(x -> x == " Note_on_c", s[:, 3])
    #next lines are to get the time in miliseconds when the pitch starts and ends.
    if nv == 1
        #nv = 1
        if isempty(n_off) || length(n_off) < length(n_on) / 2
            b = s[s[:, 1].==1, :] #takes the events of the voice i
            mat = b[b[:, 3].==" Note_on_c", :]
            mat = mat[findall(x -> x != "", mat[:, 5]), :]
            voces = get_onon_notes(mat) #construct an array of of information of initial time, finish time, pitch and intensity.
        else
            mat = s[s[:, 1].==1, :]
            mat = mat[findall(x -> x != "", mat[:, 5]), :]
            voces = get_onoff_notes(mat)
        end
    else
        voces = Array{Matrix}(undef, nv) #initialze an array
        if isempty(n_off) || length(n_off) < length(n_on) / 2
            #checks if the midi has events of note_off
            for i = 1:(nv) #if does not, it construct the series in this way
                b = s[s[:, 1].==i, :] #takes the events of the voice i
                #if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end #checks if there are notes in the channel
                mat = b[b[:, 3].==" Note_on_c", :]
                mat = mat[findall(x -> x != "", mat[:, 5]), :]
                voces[i] = get_onon_notes(mat)#construct an array of information of initial time, finish time, pitch and intensity.
            end
        else
            for i = 1:nv #if has note_off events , it does this way
                ini = findfirst(x -> x == " Note_on_c", s[s[:, 1].==i, 3]) #initial time
                fin = findlast(x -> x == " Note_off_c", s[s[:, 1].==i, 3]) #finish time
                if ini == 0 || fin == 0
                    continue
                end
                mat = s[s[:, 1].==i, :]
                mat = mat[findall(x -> x != "", mat[:, 5]), :]
                voces[i] = get_onoff_notes(mat) #construct an array of information of initial time, finish time, pitch and intensity
            end
        end
    end
    if nv != 1
        n_vs = filter_undef!(voces)
        filter!(x -> length(x) > 0, n_vs)
    else
        n_vs = []
        push!(n_vs, voces)
    end
    nv = size(n_vs)[1]
    t_max = max_tempo(n_vs, nv)
    number_measures = length(measure)
    all_v = vcat(n_vs...)
    norder = sortperm(all_v[:, 1])
    all_v = all_v[norder, :]
    all_divis = []
    if number_measures > 1
        chunked_voice = []
        for ct = 1:number_measures
            th = (q_t / qdiv) / notes_measure[ct][1]
            chunked_meas = []
            if ct < number_measures
                t_final = meas_time[ct+1] - meas_time[ct]
            else
                t_final = t_max - meas_time[ct]
            end
            n_chunks, w_res = divrem(t_final, w_s[ct])
            for i = 1:n_chunks
                st = meas_time[ct] + (i - 1) * w_s[ct]
                en = st + w_s[ct] - 1
                in_notes = findall(x -> st <= x < en, all_v[:, 1])
                st_times = all_v[in_notes, 1]
                in_times = all_v[in_notes, 1] .- st
                loc = zeros(length(in_notes))
                num_n = [sym_measure[ct] for i = 1:length(in_notes)]
                dur = all_v[in_notes, 2] - all_v[in_notes, 1]
                pitc = all_v[in_notes, 3]
                for z = 1:length(in_notes)
                    divis = in_times[z] / notes_measure[ct][1]
                    dif, bt = modf(divis)
                    if dif <= th #check if the note starts in the tempo is found, within a threshold
                        loc[z] = bt + 1
                    elseif 1 - dif <= th #if the note starts in the next tempo
                        loc[z] = bt + 2
                    elseif dif >= 0.5 - th && dif <= 0.5 + th  #check if the note starts on a half beat
                        loc[z] = bt + 1.5
                    else
                        loc[z] = 0 #the note doesn't start in the tempo, is in between.
                    end
                    push!(all_divis, dif * smdiv)
                    #println("Pitch: $(pitc[z])",'\t',"Starting time in piece (midi clock): $(st_times[z]), in the measure:$(in_times[z]) ",'\t',"Ratio Starting/LengthBeat: $(divis)",'\t',"Selected starting beat (in measure): $(loc[z])")
                end
                if abs_time
                    push!(chunked_meas, [loc num_n all_v[in_notes, 1:2] dur all_v[in_notes, 3]])
                else
                    push!(chunked_meas, [loc num_n all_v[in_notes, 3] dur])
                end
            end
            #next lines are for the last measure or the last fraction left.

            if ct == number_measures && w_res > 0
                st = meas_time[ct] + n_chunks * w_s[ct]
                en = t_max
                in_notes = findall(x -> st <= x < en, all_v[:, 1])
                st_times = all_v[in_notes, 1]
                in_times = all_v[in_notes, 1] .- st
                loc = zeros(length(in_notes))
                num_n = [sym_measure for i = 1:length(in_notes)]
                dur = all_v[in_notes, 2] - all_v[in_notes, 1]
                pitc = all_v[in_notes, 3]
                for z = 1:length(in_notes)
                    #pl = findfirst(x-> x== in_times[z], sort([in_times[z];notes_measure[ct]]))
                    divis = in_times[z] / notes_measure[1][1]
                    dif, bt = modf(divis) #- (notes_measure[ct][pl] - notes_measure[ct][1])
                    #hbeat = notes_measure[ct][1]/2
                    if dif <= th #check if the note starts in the tempo is found, within a threshold
                        loc[z] = bt + 1
                    elseif 1 - dif <= th #if the note starts in the next tempo
                        loc[z] = bt + 2
                    elseif dif >= 0.5 - th && dif <= 0.5 + th  #check if the note starts on a half beat
                        loc[z] = bt + 1.5
                    else
                        loc[z] = 0 #the note doesn't start in the tempo, is in between.
                    end
                    push!(all_divis, dif * smdiv)
                    #println("Pitch: $(pitc[z])",'\t',"Starting time in piece (midi clock): $(st_times[z]), in the measure:$(in_times[z]) ",'\t',"Ratio Starting/LengthBeat: $(divis)",'\t', "Ratio Decimals multiplied by 32: $(dif*smdiv)")#"Selected starting beat (in measure): $(loc[z])")
                end
                if abs_time
                    push!(chunked_meas, [loc num_n all_v[in_notes, 1:2] dur all_v[in_notes, 3]]) #returns also the values for the absolute time in ms when the note started and when it ended
                else
                    push!(chunked_meas, [loc num_n all_v[in_notes, 3] dur])
                end
            end
            push!(chunked_voice, chunked_meas)
        end
    else
        n_chunks, w_res = divrem(t_max, w_s)
        th = (q_t / qdiv) / notes_measure[1][1]
        chunked_voice = []
        for i = 1:n_chunks
            st = (i - 1) * w_s
            en = st + w_s - 1
            in_notes = findall(x -> st <= x < en, all_v[:, 1])
            st_times = all_v[in_notes, 1]
            in_times = all_v[in_notes, 1] .- st
            loc = zeros(length(in_notes))
            num_n = [sym_measure for i = 1:length(in_notes)]
            dur = all_v[in_notes, 2] - all_v[in_notes, 1]
            pitc = all_v[in_notes, 3]
            for z = 1:length(in_notes)
                divis = in_times[z] / notes_measure[1][1]
                dif, bt = modf(divis)
                if dif <= th #check if the note starts in the tempo is found, within a threshold
                    loc[z] = bt + 1
                elseif 1 - dif <= th #if the note starts in the next tempo
                    loc[z] = bt + 2
                elseif dif >= 0.5 - th && dif <= 0.5 + th  #check if the note starts on a half beat
                    loc[z] = bt + 1.5
                else
                    loc[z] = 0 #the note doesn't start in the tempo, is in between.
                end
                push!(all_divis, dif * smdiv)
                #println("Pitch: $(pitc[z])",'\t',"Starting time in piece (midi clock): $(st_times[z]), in the measure:$(in_times[z]) ",'\t',"Ratio Starting/LengthBeat: $(divis)",'\t',"Selected starting beat (in measure): $(loc[z])")
            end
            if abs_time
                push!(chunked_voice, [loc num_n all_v[in_notes, 1:2] dur all_v[in_notes, 3]])
            else
                push!(chunked_voice, [loc num_n all_v[in_notes, 3] dur])
            end

        end
        if w_res > 0
            st = n_chunks * w_s
            if st != t_max
                en = t_max
                #println(st,"\t",en)
                in_notes = findall(x -> st <= x < en, all_v[:, 1])
                st_times = all_v[in_notes, 1]
                in_times = all_v[in_notes, 1] .- st
                loc = zeros(length(in_notes))
                num_n = [sym_measure for i = 1:length(in_notes)]
                dur = all_v[in_notes, 2] - all_v[in_notes, 1]
                pitc = all_v[in_notes, 3]
                for z = 1:length(in_notes)
                    #pl = findfirst(x-> x== in_times[z], sort([in_times[z];notes_measure[ct]]))
                    divis = in_times[z] / notes_measure[1][1]
                    dif, bt = modf(divis) #- (notes_measure[ct][pl] - notes_measure[ct][1])
                    #hbeat = notes_measure[ct][1]/2
                    if dif <= th #check if the note starts in the tempo is found, within a threshold
                        loc[z] = bt + 1
                    elseif 1 - dif <= th #if the note starts in the next tempo
                        loc[z] = bt + 2
                    elseif dif >= 0.5 - th && dif <= 0.5 + th  #check if the note starts on a half beat
                        loc[z] = bt + 1.5
                    else
                        loc[z] = 0 #the note doesn't start in the tempo, is in between.
                    end
                    push!(all_divis, dif * smdiv)
                    #println("Pitch: $(pitc[z])",'\t',"Starting time in piece (midi clock): $(st_times[z]), in the measure:$(in_times[z]) ",'\t',"Ratio Starting/LengthBeat: $(divis)",'\t', "Ratio Decimals multiplied by 32: $(dif*smdiv)")#"Selected starting beat (in measure): $(loc[z])")
                end
                if abs_time
                    push!(chunked_voice, [loc num_n all_v[in_notes, 1:2] dur all_v[in_notes, 3]])
                else
                    push!(chunked_voice, [loc num_n all_v[in_notes, 3] dur])
                end
            end
        end
    end
    #next lines are for the last measure or the last fraction left.
    residues = map(x -> modf(x)[1], all_divis)
    n_off = length(findall(x -> x > 10^-6, residues))
    res_frac = n_off / length(residues)
    if number_measures > 1
        out_measures = vcat(chunked_voice...)
    else
        out_measures = chunked_voice
    end
    #if res_frac > 0.05
    #   println("WARNING!!! \n THE FRACTION OF NOTES FALLING OUTSIDE THE MIDI CLOCK IS: $(res_frac)")
    #end
    return filter(x -> !isempty(x), out_measures), [length(residues) n_off]
end

"""
    function get_piece_by_measure_m21(m21_data; start_measure=1)

    Returns an array with the information of pitches and duration on each measure in the following format:

    | #measure | time signature | start_quarter | end_quarter | duration (in quarters) | pitch (0-127) |
"""
function get_piece_by_measure_m21(m21_data; start_measure=1)

    notes = []
    for parts in m21_data.parts #getting all the notes in all parts, with duration and pitch
        for note in parts.flat.notes
            try
                local start, duration
                measure = note.measureNumber + (start_measure - 1)
                start = float(note.offset)
                duration = float(note.quarterLength)
                pitch = note.pitch.ps
                #print(start,'\n')
                push!(notes, [measure start start + duration duration pitch])
            catch
                local start, duration
                measure = note.measureNumber + (start_measure - 1)
                try

                    start = float(note.offset.())
                catch
                    try
                        numer = note.offset.numerator
                        denom = note.offset.denominator
                        start = float(numer / denom)
                    catch
                        start = float(note.offset)
                    end
                end
                #print(note.quarterLength,'\n')
                try
                    #print("here 1 \n")
                    duration = float(note.quarterLength.())
                catch
                    try
                        #print("here 2 \n")
                        numer = note.quarterLength.numerator
                        denom = note.quarterLength.denominator
                        duration = float(numer / denom)
                    catch
                        #print("here 3 \n")
                        duration = float(note.quarterLength)

                    end
                end
                #print("duration before appending: ",duration, '\n')
                #print("start before appending: ", start, '\n')
                for chord_note in note.pitches
                    pitch = chord_note.ps
                    #print("aqui \n")
                    #print(start,'\t',duration,'\n')
                    push!(notes, [measure start start + duration duration pitch])

                end
            end
        end
    end
    all_notes = vcat(notes...)
    n_measures = Int(maximum(all_notes[:, 1]))
    timesig = []
    lquarters = []
    for i in 0:n_measures #extracting time signatures
        try
            ts = m21_data.parts[1].measure(i).timeSignature
            push!(timesig, [i ts.ratioString])
            nquarters = ts.numerator / ts.denominator * 4
            push!(lquarters, [i nquarters])
        catch
        end
    end

    measures_piece = []
    if length(timesig) > 1
        for n_ts in 2:length(timesig)
            last_m = timesig[n_ts][end, 1] #get the last measure where the time signature applies 
            first_m = timesig[n_ts-1][1, 1]
            ix = findall(x -> first_m <= x <= last_m, all_notes[:, 1])
            m_ts = [timesig[n_ts-1][1, 2] for i in 1:length(ix)]
            push!(measures_piece, [all_notes[ix, 1] m_ts all_notes[ix, 2:end]])
        end
        #now last time signature
        last_m = n_measures
        first_m = timesig[end][1, 1]
        ix = findall(x -> first_m <= x <= last_m, all_notes[:, 1])
        m_ts = [timesig[end][1, 2] for i in 1:length(ix)]
        push!(measures_piece, [all_notes[ix, 1] m_ts all_notes[ix, 2:end]])
    else
        ix = all_notes[:, 1]
        m_ts = [timesig[1][1, 2] for i in 1:length(ix)]
        push!(measures_piece, [all_notes[:, 1] m_ts all_notes[:, 2:end]])
    end
    out_measures = vcat(measures_piece...)
    return out_measures[sortperm(out_measures[:, 1]), :]
end

"""
    get_piece_by_measure_mid(filepath; abs_time=true, start_measure=1, ignore_drums=true)

    Parse a raw `.mid` file with MIDI.jl and return a flat `N×6` matrix in the same
    layout as `get_piece_by_measure_m21`:

    | #measure | time signature | start | end | duration | pitch (0-127) |

    The time base is tick-exact (derived from the file's ticks-per-quarter and its
    time-signature / tempo meta events), so no `qdiv` beat-snapping tolerance is needed.

    - `abs_time=true`  → start/end/duration are in **milliseconds** (matches the CSV path,
      honoring all tempo changes).
    - `abs_time=false` → start/end/duration are in **absolute quarter notes** (matches the
      XML path); this is tempo-independent.

    `start_measure` offsets the reported measure numbers (matching `get_piece_by_measure_m21`).
    `ignore_drums=true` drops MIDI channel 10 (0-indexed channel 9), whose pitches are
    percussion-map indices rather than real pitches.
"""
function get_piece_by_measure_mid(filepath; abs_time=true, start_measure=1, ignore_drums=true)
    midi = MIDI.load(filepath)
    tpq = Int(midi.tpq)                       # ticks per quarter note

    # --- 1. Collect time-signature and tempo meta events (absolute ticks) ---
    tsigs = Tuple{Int,Int,Int}[]              # (tick, numerator, denominator)
    tempos = Tuple{Int,Int}[]                 # (tick, microseconds_per_quarter)
    for tr in midi.tracks
        t = 0
        for ev in tr.events
            t += ev.dT
            if ev isa MIDI.TimeSignatureEvent
                push!(tsigs, (t, ev.numerator, ev.denominator))
            elseif ev isa MIDI.SetTempoEvent
                push!(tempos, (t, ev.tempo))
            end
        end
    end
    isempty(tsigs) && push!(tsigs, (0, 4, 4))             # default 4/4
    isempty(tempos) && push!(tempos, (0, 500_000))        # default 120 qpm
    sort!(tsigs); sort!(tempos)

    # --- 2. Gather notes from every track (handles note_off and note_on(vel=0)) ---
    all_notes = MIDI.Note[]
    for tr in midi.tracks
        append!(all_notes, MIDI.getnotes(tr).notes)
    end
    ignore_drums && filter!(n -> n.channel != 9, all_notes)
    sort!(all_notes; by = n -> n.position)
    isempty(all_notes) && return Matrix{Any}(undef, 0, 6)

    # --- 3. Build measure boundaries from time-signature events + tpq ---
    #   ticks_per_measure = tpq * 4 * num / den
    end_tick = maximum(n.position + n.duration for n in all_notes)
    measure_bounds = Tuple{Int,Int,Int,Int}[]   # (measure_no, start_tick, end_tick, ts_index)
    measure_no = 1
    cursor = 0
    for (i, (ts_tick, num, den)) in enumerate(tsigs)
        next_ts_tick = i < length(tsigs) ? tsigs[i+1][1] : end_tick + 1
        mlen = div(tpq * 4 * num, den)
        while cursor < next_ts_tick && cursor <= end_tick
            push!(measure_bounds, (measure_no, cursor, cursor + mlen, i))
            measure_no += 1
            cursor += mlen
        end
    end

    # --- 4. Tick → ms helper (piecewise-constant tempo) ---
    function tick_to_ms(tick)
        ms = 0.0
        for (i, (tt, μspq)) in enumerate(tempos)
            next_tt = i < length(tempos) ? tempos[i+1][1] : tick
            seg_end = min(tick, next_tt)
            ms += (seg_end - tt) * μspq / tpq / 1000
            tick <= next_tt && break
        end
        return ms
    end

    # --- 5. Assign each note to a measure and build the output matrix ---
    rows = Vector{Any}[]
    for n in all_notes
        idx = findfirst(b -> b[2] <= n.position < b[3], measure_bounds)
        idx === nothing && (idx = length(measure_bounds))   # clamp to last measure
        m_no = measure_bounds[idx][1] + (start_measure - 1)
        ts_i = measure_bounds[idx][4]
        ts_str = "$(tsigs[ts_i][2])/$(tsigs[ts_i][3])"
        if abs_time
            s = tick_to_ms(n.position)
            e = tick_to_ms(n.position + n.duration)
        else
            s = n.position / tpq
            e = (n.position + n.duration) / tpq
        end
        push!(rows, Any[m_no, ts_str, s, e, e - s, Int(n.pitch)])
    end
    return permutedims(reduce(hcat, rows))       # N × 6 matrix
end

const STEP_TO_SEMITONE = Dict(
    "C" => 0, "D" => 2, "E" => 4, "F" => 5, "G" => 7, "A" => 9, "B" => 11
)

"""
    pitch_to_midi(pitch_elem) -> Int

Convert a MusicXML `<pitch>` element to a MIDI note number (matching music21's `.ps`):

    MIDI = 12 * (octave + 1) + step_semitone + alter
"""
function pitch_to_midi(pitch_elem)
    step = nodecontent(findfirst("step", pitch_elem))
    octave = parse(Int, nodecontent(findfirst("octave", pitch_elem)))
    alter_node = findfirst("alter", pitch_elem)
    alter = alter_node === nothing ? 0 : round(Int, parse(Float64, nodecontent(alter_node)))
    return 12 * (octave + 1) + STEP_TO_SEMITONE[step] + alter
end

"""
    read_musicxml(filepath) -> EzXML.Document

Read a `.musicxml`, `.xml`, or `.mxl` file and return the parsed XML document.
`.mxl` files are ZIP archives; the inner score is located by reading the
`full-path` of `<rootfile>` in `META-INF/container.xml` (with a fallback to the
first non-`META-INF` XML entry).
"""
function read_musicxml(filepath)
    if endswith(lowercase(filepath), ".mxl")
        zr = ZipFile.Reader(filepath)
        try
            fullpath = nothing
            cidx = findfirst(f -> f.name == "META-INF/container.xml", zr.files)
            if cidx !== nothing
                cdoc = parsexml(read(zr.files[cidx], String))
                rf = findfirst("//rootfile", root(cdoc))
                rf !== nothing && (fullpath = rf["full-path"])
            end
            if fullpath === nothing
                idx = findfirst(f -> !startswith(f.name, "META-INF") &&
                        (endswith(lowercase(f.name), ".xml") ||
                         endswith(lowercase(f.name), ".musicxml")), zr.files)
                idx === nothing && error("No MusicXML entry found inside $filepath")
                fullpath = zr.files[idx].name
            end
            fidx = findfirst(f -> f.name == fullpath, zr.files)
            fidx === nothing && error("Rootfile $fullpath not found inside $filepath")
            return parsexml(read(zr.files[fidx], String))
        finally
            close(zr)
        end
    else
        return readxml(filepath)
    end
end

"""
    get_piece_by_measure_xml(filepath; start_measure=1)

Parse a MusicXML file (`.musicxml`, `.xml`, or `.mxl`) natively with EzXML and
return a flat `N×6` matrix in the same layout as `get_piece_by_measure_m21`:

    | #measure | time signature | start_quarter | end_quarter | duration (quarters) | pitch (0-127) |

Notes from all parts are merged; the time-signature timeline is taken from the
first part (matching the previous music21 behavior). Durations are exact in
quarter notes (`<duration>/<divisions>`); chords share an onset, rests and
unpitched notes are skipped, tied noteheads are kept separate, and `<forward>`/
`<backup>` shift the within-measure offset for multi-voice parts.
"""
function get_piece_by_measure_xml(filepath; start_measure=1)
    doc = read_musicxml(filepath)
    r = root(doc)
    nodename(r) == "score-timewise" &&
        error("score-timewise MusicXML is not supported; convert to score-partwise.")
    parts = findall("//part", r)
    isempty(parts) && return Matrix{Any}(undef, 0, 6)

    # --- Time-signature timeline from the first part (measure_number => "n/d") ---
    ts_changes = Tuple{Int,String}[]
    for measure in findall("measure", parts[1])
        tn = findfirst(".//time", measure)
        if tn !== nothing
            bnode = findfirst("beats", tn)
            btnode = findfirst("beat-type", tn)
            if bnode !== nothing && btnode !== nothing
                push!(ts_changes,
                    (parse(Int, measure["number"]), nodecontent(bnode) * "/" * nodecontent(btnode)))
            end
        end
    end
    function ts_for(m)
        ts = "4/4"
        for (mn, s) in ts_changes
            mn <= m ? (ts = s) : break
        end
        return ts
    end

    # --- Walk every part, accumulating absolute quarter offsets ---
    rows = Vector{Any}[]
    for part in parts
        divisions = 1
        measure_base = 0.0          # absolute quarter offset at the start of this measure
        for measure in findall("measure", part)
            raw_mnum = parse(Int, measure["number"])
            measure_num = raw_mnum + (start_measure - 1)
            ts_str = ts_for(raw_mnum)
            local_offset = 0.0
            max_offset = 0.0
            prev_onset = 0.0
            for child in eachelement(measure)
                tag = nodename(child)
                if tag == "attributes"
                    dn = findfirst("divisions", child)
                    dn !== nothing && (divisions = parse(Int, nodecontent(dn)))
                elseif tag == "forward"
                    dn = findfirst("duration", child)
                    dn !== nothing && (local_offset += parse(Float64, nodecontent(dn)) / divisions)
                    max_offset = max(max_offset, local_offset)
                elseif tag == "backup"
                    dn = findfirst("duration", child)
                    dn !== nothing && (local_offset -= parse(Float64, nodecontent(dn)) / divisions)
                elseif tag == "note"
                    dn = findfirst("duration", child)
                    dur_q = dn === nothing ? 0.0 : parse(Float64, nodecontent(dn)) / divisions
                    if findfirst("chord", child) !== nothing
                        onset = prev_onset                  # chord notes share the previous onset
                    else
                        onset = local_offset
                        prev_onset = onset
                        local_offset += dur_q
                        max_offset = max(max_offset, local_offset)
                    end
                    pitch_elem = findfirst("pitch", child)
                    (findfirst("rest", child) !== nothing || pitch_elem === nothing) && continue
                    abs_onset = measure_base + onset
                    push!(rows, Any[measure_num, ts_str, abs_onset, abs_onset + dur_q,
                        dur_q, pitch_to_midi(pitch_elem)])
                end
            end
            measure_base += max_offset
        end
    end
    isempty(rows) && return Matrix{Any}(undef, 0, 6)
    out = permutedims(reduce(hcat, rows))
    return out[sortperm(out[:, 1]), :]                      # sort by measure, matching m21
end

"""
    function get_piece_by_measure(data; csv=true, qdiv=32, abs_time=true, start_measure=1)

    Dispatch to the right parser based on the input. A `.mid` path uses the raw-MIDI
    parser; a `.musicxml` / `.xml` / `.mxl` path uses the native EzXML parser, falling
    back to music21 (via PyCall) if the native parser errors or returns no notes;
    otherwise `csv=true` selects the midicsv parser and `csv=false` the (legacy)
    music21 parser for a pre-parsed music21 `Score` object.
"""

function get_piece_by_measure(data; csv=true, qdiv=32, abs_time=true, start_measure=1)
    if data isa AbstractString && endswith(lowercase(data), ".mid")
        return get_piece_by_measure_mid(data, abs_time=abs_time, start_measure=start_measure)
    elseif data isa AbstractString && (endswith(lowercase(data), ".musicxml") ||
                                       endswith(lowercase(data), ".mxl") ||
                                       endswith(lowercase(data), ".xml"))
        # Native Julia parser first; fall back to music21 if it errors or finds nothing.
        local piece
        try
            piece = get_piece_by_measure_xml(data, start_measure=start_measure)
        catch e
            @warn "Native XML parser failed; falling back to music21." path=data exception=e
            return _get_piece_by_measure_m21_from_path(data, start_measure=start_measure)
        end
        if size(piece, 1) == 0
            @warn "Native XML parser returned no notes; falling back to music21." path=data
            return _get_piece_by_measure_m21_from_path(data, start_measure=start_measure)
        end
        return piece
    elseif csv
        return get_piece_by_measure_csv(data, qdiv=qdiv, abs_time=abs_time)
    else
        return get_piece_by_measure_m21(data, start_measure=start_measure)
    end
end

"""
    _get_piece_by_measure_m21_from_path(path; start_measure=1)

Fallback helper: parse a score file with music21 (via PyCall) and run the legacy
`get_piece_by_measure_m21`. Used only when the native EzXML parser fails.
"""
function _get_piece_by_measure_m21_from_path(path; start_measure=1)
    m21 = pyimport("music21")
    score = m21.converter.parse(path)
    return get_piece_by_measure_m21(score, start_measure=start_measure)
end

"""
    get_center_effect(chunk_notes; r=1, h=sqrt(2/15), 
                    all_keys=all_keys, pos_all_keys=pos_all_keys, 
                    sbeat_w=[[1.],[1.]], lin_w=1)
    
    Computes the coordinates (x,y,z) for the center of effect (mass)
    for a given set of notes, these notes can be in different formats:
    Notes are always represented in MIDI notation, note durations can
    be in any numeric format.

        Matrix: From the DataFrame built with get_xml_df or get_csv_df
        1d - Array: [n1,n2,...,nm], where each element represents a note in MIDI notation
        two 1d - Array: [n1,n2,...,nm], [d1,d2,...,dm], with the set of n
        representing the notes in MIDI notation and d each of their durations.

"""
function get_center_effect(chunk_notes::Matrix{<:Any}; r=1, h=sqrt(2 / 15), all_keys=all_keys, pos_all_keys=pos_all_keys, sbeat_w=[[1.0], [1.0]], lin_w=1)
    ptcs = convert(Vector{Int},chunk_notes[:, 6])
    durs = chunk_notes[:, 5]
    pbeat = chunk_notes[:, 1]
    beat_w = ones(length(durs)) #array of the beat weights
    for b = 1:length(sbeat_w[1])
        loc_b = findall(x -> x == sbeat_w[1][b], pbeat) #finding all notes that start at beat sbeat_w[1][b]
        beat_w[loc_b] .= sbeat_w[2][b] #this is the weight.
    end
    notas, n_we = get_local_lin_w(ptcs, lin_w) #doing the linear weight in the pitches
    ii = vcat(map(x -> findall(y -> y == x, notas), ptcs)...)
    #println(ptcs)
    b_wei = n_we[ii] #getting the linear weight for every note i n the array of pitches

    #spi_ix = cluster_notes(ptcs) .+ 24
    #map to spiral array
    spi_notes = map(x -> get_cfpitch(x), ptcs)
    spi_ix = circular_to_linear(spi_notes)
    spi_p = map(x -> get_pitch(x, r=r, h=h), spi_ix) #getting the location (x,y,z) for each pitch
    t_ws = map((x, y, z) -> x * y * z, beat_w, durs, b_wei) #computing the total weights
    cv_i = map((x, y) -> x * y, t_ws, spi_p) / sum(t_ws) #computing the location of the pitches with their relative weights
    c_i = sum(cv_i) #finding the center of effect
    return c_i
end

function get_center_effect(seq_notes::AbstractVector{<:Number}; r=1, h=sqrt(2 / 15), all_keys=all_keys, pos_all_keys=pos_all_keys, sbeat_w=[[1.0], [1.0]], lin_w=1)
    ptcs = convert(Vector{Int},seq_notes)
    durs = [1 for _ = 1:length(seq_notes)]
    pbeat = [1 for _ = 1:length(seq_notes)]
    beat_w = ones(length(durs)) #array of the beat weights
    for b = 1:length(sbeat_w[1])
        loc_b = findall(x -> x == sbeat_w[1][b], pbeat) #finding all notes that start at beat sbeat_w[1][b]
        beat_w[loc_b] .= sbeat_w[2][b] #this is the weight.
    end
    notas, n_we = get_local_lin_w(ptcs, lin_w) #doing the linear weight in the pitches
    ii = vcat(map(x -> findall(y -> y == x, notas), ptcs)...)
    #println(ptcs)
    b_wei = n_we[ii] #getting the linear weight for every note i n the array of pitches

    #spi_ix = cluster_notes(ptcs) .+ 24
    spi_notes = map(x -> get_cfpitch(x), ptcs)
    spi_ix = circular_to_linear(spi_notes)
    spi_p = map(x -> get_pitch(x, r=r, h=h), spi_ix) #getting the location (x,y,z) for each pitch
    t_ws = map((x, y, z) -> x * y * z, beat_w, durs, b_wei) #computing the total weights
    cv_i = map((x, y) -> x * y, t_ws, spi_p) / sum(t_ws) #computing the location of the pitches with their relative weights
    c_i = sum(cv_i) #finding the center of effect
    return c_i
end

function get_center_effect(notes::AbstractVector{<:Number}, durs::AbstractVector{<:Number}; r=1, h=sqrt(2 / 15), all_keys=all_keys, pos_all_keys=pos_all_keys, sbeat_w=[[1.0], [1.0]], lin_w=1)
    ptcs = Int.(notes)
    #durs = [1 for i = 1:length(seq_notes)]
    pbeat = [1 for _ = 1:length(notes)]
    beat_w = ones(length(durs)) #array of the beat weights
    for b = 1:length(sbeat_w[1])
        loc_b = findall(x -> x == sbeat_w[1][b], pbeat) #finding all notes that start at beat sbeat_w[1][b]
        beat_w[loc_b] .= sbeat_w[2][b] #this is the weight.
    end
    notas, n_we = get_local_lin_w(ptcs, lin_w) #doing the linear weight in the pitches
    ii = vcat(map(x -> findall(y -> y == x, notas), ptcs)...)
    #println(ptcs)
    b_wei = n_we[ii] #getting the linear weight for every note i n the array of pitches

    #spi_ix = cluster_notes(ptcs) .+ 24
    spi_notes = map(x -> get_cfpitch(x), ptcs)
    spi_ix = circular_to_linear(spi_notes)
    
    spi_p = map(x -> get_pitch(x, r=r, h=h), spi_ix) #getting the location (x,y,z) for each pitch
    t_ws = map((x, y, z) -> x * y * z, beat_w, durs, b_wei) #computing the total weights
    cv_i = map((x, y) -> x * y, t_ws, spi_p) / sum(t_ws) #computing the location of the pitches with their relative weights
    c_i = sum(cv_i) #finding the center of effect
    return c_i
end

function get_distance_to_keys(c_i)
    d_to_keys = round.(map(x -> euclidean(c_i, x), pos_all_keys), digits=4) #computing the eclidean distance to all keys

    ranking = sortperm(d_to_keys) #ranking the distances from the closest to the farthest

    return [all_keys[ranking] d_to_keys[ranking]][1:22, :]
end

function get_distance_to_chords(c_i)
    d_to_chords = round.(map(x -> euclidean(c_i, x), pos_all_chords), digits=4)
    ranking = sortperm(d_to_chords)
    return [all_chords[ranking] d_to_chords[ranking]][1:22, :]
end

function get_xml_df(piece_xml)
    piece = get_piece_by_measure(piece_xml, csv=false)
    df_piece = DataFrame(
        :Measure => convert(Array{Int64,1}, piece[:, 1]),
        :TimeSignature => piece[:, 2],
        :StartQuarter => piece[:, 3],
        :EndQuarter => piece[:, 4],
        :Duration => piece[:, 5],
        :Pitch => convert(Array{Int64,1}, piece[:, 6])
    )
    return df_piece
end
function get_mid_df(filepath; abs_time=true, start_measure=1, ignore_drums=true)
    piece = get_piece_by_measure_mid(filepath, abs_time=abs_time,
        start_measure=start_measure, ignore_drums=ignore_drums)
    df_piece = DataFrame(
        :Measure => convert(Array{Int64,1}, piece[:, 1]),
        :TimeSignature => piece[:, 2],
        (abs_time ? :StartTime : :StartQuarter) => piece[:, 3],
        (abs_time ? :EndTime : :EndQuarter) => piece[:, 4],
        :Duration => piece[:, 5],
        :Pitch => convert(Array{Int64,1}, piece[:, 6])
    )
    return df_piece
end
function get_csv_df(piece_csv)
    piece = get_piece_by_measure(piece_csv, csv=true)[1]
    num_mea = []
    for nm in 1:length(piece)
        push!(num_mea, [nm for i in 1:size(piece[nm], 1)])
    end
    num_mea = vcat(num_mea...)
    piece = vcat(piece...)
    df_piece = DataFrame(
        :Measure => convert(Array{Int64,1}, num_mea),
        :TimeSignature => piece[:, 2],
        :StartTime => piece[:, 3],
        :EndTime => piece[:, 4],
        :Duration => piece[:, 5],
        :Pitch => convert(Array{Int64,1}, piece[:, 6])
    )
    return df_piece
end