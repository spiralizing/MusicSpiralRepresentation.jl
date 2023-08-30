global pitch_names = ["C","G","D","A","E","B","Gb/F#","Db","Ab","Eb","Bb","F"]
global nMajor_keys = ["C", "G", "D", "A","E", "B/Cb","Gb/F#","Db/C#","Ab","Eb", "Bb","F"]
global nminor_keys = ["c","g","d","a","e","b","f#","c#","g#/ab","eb/d#","bb","f"]
#handwritting the functional harmony notations "key or chord relative to the fundamental"
#Major roman numerals, for the coe notation.
Major_RN = Dict(0 => "I",
    7 => "#I/bII",
    2 => "II",
    9 => "#II/bIII",
    4 => "III",
    #"III#/IVb" = 11
    11 => "IV",
    6 => "#IV/bV",
    1 => "V",
    8 => "#V/bVI",
    3 => "VI",
    10 => "#VI/bVII",
    5 => "VII",
)
##--
Minor_RN = Dict(0 => "i",
    7 => "#i/bii",
    2 => "ii",
    9 => "#ii/biii",
    4 => "iii",
    #"III#/IVb" = 11
    11 => "iv",
    6 => "#iv/bv",
    1 => "v",
    8 => "#v/bvi",
    3 => "vi",
    10 => "#vi/bvii",
    5 => "vii",
)
global nall_keys = [] #All keys
for (a,b) in zip(nminor_keys,nMajor_keys)
    push!(nall_keys,a); push!(nall_keys,b)
end
global cf_notes = [0, 7, 2, 9, 4, 11, 6, 1, 8, 3, 10, 5]
global midi_notes = ["C","C#","D","Eb","E","F","F#","G","G#","A","Bb","B"]
#the real value of the notes in the circle of fifths, chromatic scale representation from C = 0 to B = 11
# r = 1
# h = sqrt(2/15)
# w = [0.536, 0.274, 0.19] #weights used for the center of effect algorithm
# u = v = o = w
# a = 0.75
# b = 0.75
global h_octav = 4.3818

"""
    get_pitch(k; r=1, h=sqrt(2/15))

    Returns the location (x,y,z) of a pitch "k"
    in the spiral representation, pitches are in module 12
    the variable "r" is the radius of the cylinder
    and "h" is the vertical distance between thirds.

"""
function get_pitch(k; r=1, h=sqrt(2 / 15))
    p = [r * sin(k * pi / 2), r * cos(k * pi / 2), k * h]
    return p
end

"""
    get_Major_chord(k; w=[0.536, 0.274, 0.19], r=1, h=sqrt(2/15))

    Returns the location (x,y,z) of a Major chord given by a linear combination
    of a triad constructed from a fundamental note "k", the fifth (k+1) and the third (k+4)
    in the spiral representation. The vector "w" contains the weights used for each term.
"""
function get_Major_chord(k; w=[0.536, 0.274, 0.19], r=1, h=sqrt(2 / 15))
    CM = w[1] * get_pitch(k, r=r, h=h) + w[2] * get_pitch(k + 1, r=r, h=h) + w[3] * get_pitch(k + 4, r=r, h=h)
    return CM
end

"""
    get_minor_chord(k; u=[0.536, 0.274, 0.19], r=1, h=sqrt(2/15))

    Returns the location (x,y,z) of a Minor chord given by a linear combination
    of a triad constructed from a fundamental note "k", the fifth (k+1) and the minor third (k-3)
    in the spiral representation. The vector "u" contains the weights used for each term.
"""
function get_minor_chord(k; u=[0.536, 0.274, 0.19], r=1, h=sqrt(2 / 15))
    Cm = u[1] * get_pitch(k, r=r, h=h) + u[2] * get_pitch(k + 1, r=r, h=h) + u[3] * get_pitch(k - 3, r=r, h=h)
    return Cm
end

"""
    get_Major_key(k; o=[0.536, 0.274, 0.19], w=[0.536, 0.274, 0.19], r=1, h=sqrt(2/15))

    Returns the location (x,y,z) of a Major Key given by a linear combination
    of chord triads constructed from the major chords of a fundamental "tonic" note "k",
    the dominant/fifth (k+1) and the subdominant/minor third (k-3) in the spiral representation.
    The vector "o" contains the weights used for each chord and "w" are the weights on each note in the chord.
"""
function get_Major_key(k; o=[0.536, 0.274, 0.19], w=[0.536, 0.274, 0.19], r=1, h=sqrt(2 / 15))
    TM = o[1] * get_Major_chord(k, w=w, r=r, h=h) + o[2] * get_Major_chord(k + 1, w=w, r=r, h=h) + o[3] * get_Major_chord(k - 1, w=w, r=r, h=h)
    return TM
end

"""
    get_minor_key(k; v=[0.536, 0.274, 0.19], w=[0.536, 0.274, 0.19],
    u=[0.536, 0.274, 0.19], a=0.75, b=0.75, r=1, h=sqrt(2/15))

    Returns the location (x,y,z) of a Minor Key given by a linear combination
    of chord triads constructed from the major chords of a fundamental note "k",
    the possible "dominant" triad (k+1) and the possible "subdominant" triad (k-3) in the spiral representation.
    The vector "v" contains the weights used for each chord and "w" are the weights on each note in the major chords.
    "u" are the weights for the minor chords. The dominant and subdominant minor/major chords are weighted by "a" and "b".
"""
function get_minor_key(k; v=[0.536, 0.274, 0.19], w=[0.536, 0.274, 0.19], u=[0.536, 0.274, 0.19], a=0.75, b=0.75, r=1, h=sqrt(2/15))
    Tm = v[1]*get_minor_chord(k, u=u, r=r, h=h) + v[2]*(a*get_Major_chord(k+1, w=w, r=r, h=h) + (1-a)*get_minor_chord(k+1, u=u, r=r, h=h)) + v[3]*(b*get_minor_chord(k-1, u=u, r=r, h=h) + (1-b)*get_Major_chord(k-1, w=w, r=r, h=h))
    return Tm
end

notes = [i for i = 0:127]
pitches = map(x -> get_pitch(x), notes)
major_chords = map(x -> get_Major_chord(x), notes)
minor_chords = map(x -> get_minor_chord(x), notes)
major_keys = map(x -> get_Major_key(x), notes)
minor_keys = map(x -> get_minor_key(x), notes)
pos_all_keys = vcat(major_keys, minor_keys)
#Defining all circle of fifth notes and names of major and minor keys
midi_note_names = []
all_cf_notes = []
all_Major_keys = []
all_minor_keys = []
all_pitch_names = []
push!(all_cf_notes, cf_notes)
push!(all_Major_keys, nMajor_keys)
push!(all_minor_keys, nminor_keys)
push!(midi_note_names, midi_notes)
push!(all_pitch_names, pitch_names)

for i = 1:9
    push!(all_cf_notes, cf_notes .+ (12 * i))
    push!(all_Major_keys, nMajor_keys)
    push!(all_minor_keys, nminor_keys)
    push!(midi_note_names, midi_notes)
    push!(all_pitch_names, pitch_names)
end
push!(all_cf_notes, cf_notes[1:8] .+ (120))
push!(all_Major_keys, nMajor_keys[1:8])
push!(all_minor_keys, nminor_keys[1:8])
push!(midi_note_names, midi_notes[1:8])
push!(all_pitch_names, pitch_names[1:8])


all_cf_notes = vcat(all_cf_notes...)
all_Major_keys = vcat(all_Major_keys...)
all_minor_keys = vcat(all_minor_keys...)
all_keys = vcat(all_Major_keys, all_minor_keys)
midi_note_names = vcat(midi_note_names...)
all_pitch_names = vcat(all_pitch_names...)

"""
    get_cfpitch_mod12(pitch_seq)

    Converts a pitch sequence of MIDI values (0-127) into a pitch sequence in mod12 ordered by fifths
    to be consistent with the cylindrical representation, starting in C=0, G=1, ..., etc.
"""
function get_cfpitch_mod12(pitch_seq)
    m12v = convert(Array{Int64,1},map(x->mod(x,12),pitch_seq))
    return map(y->findfirst(x->x==y,cf_notes)-1,m12v)
end

"""
    get_cfpitch(pitch_seq)

    Converts a pitch sequence of MIDI values (0-127) into a pitch sequence ordered by fifths
    to be consistent with the cylindrical representation, starting in C=0, G=1, ..., etc.
"""
function get_cfpitch(pitch_seq)
    return map(y->findfirst(x->x==y,all_cf_notes)-1,pitch_seq)
end