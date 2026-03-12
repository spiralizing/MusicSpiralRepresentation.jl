function get_onoff_notes(s)
    ini = s[s[:, 3].==" Note_on_c", :] #getting the the starting places
    fin = s[s[:, 3].==" Note_off_c", :] #ending
    dn = unique(ini[:, 5]) #different notes.
    ndn = length(dn) #number of different notes
    posin_dn = Array{Vector}(undef, ndn) #initialize arrays for getting information of each note
    posfn_dn = Array{Vector}(undef, ndn)
    in_fi = Array{Matrix}(undef, ndn)
    for i = 1:ndn       #for each different note it gets when are they played
        tmp = ini[ini[:, 5].==dn[i], :][:, 2]
        tmp2 = fin[fin[:, 5].==dn[i], :][:, 2]
        l1 = length(tmp)
        l2 = length(tmp2)
        if l1 != l2
            nt = min(l1, l2)
            note = [dn[i] for j = 1:nt]
            in_fi[i] = [tmp[1:nt] tmp2[1:nt] note]
        else
            note = [dn[i] for j = 1:length(tmp)]
            in_fi[i] = [tmp tmp2 note]
        end
    end
    nmat = vcat(in_fi...) #getting all notes in the same array
    return nmat[sortperm(nmat[:, 1]), :]  #returns the array sorted by appereance.
end

function get_new_voice(v)
    new_v = []
    ids = []
    iids = [i for i = 1:size(v)[1]]
    c = 1
    while c <= size(v)[1] - 1
        if v[c+1, 1] <= v[c, 1]
            #println(v[c,:],'\t',v[c+1,:])
            push!(new_v, v[c+1, :])
            push!(ids, c + 1)
            c += 1
        else
            #push!(ids,c)
            c += 1
        end
        #c += 1
        #println(c)
    end
    deleteat!(iids, ids)
    old_v = v[iids, :]
    new_v = convert(Array{Any,2}, transpose(reshape(vcat(new_v...), 3, :)))
    return [old_v, new_v]
end

function find_more_voices(v)
    v_1 = []
    nvs = get_new_voice(v)
    push!(v_1, nvs[1])
    lv = nvs[2]
    br = false
    while br == false #this is to check if there are more voices
        ns = get_new_voice(lv)
        push!(v_1, ns[1])
        br = isempty(ns[2])
        if br
            break
        end
        #println(br)
        lv = ns[2]
    end
    return v_1
end

function filter_undef!(voces)
    tt = map(x -> isassigned(x), voces)
    deleteat!(voces, findall(x -> x == false, tt))
end
function max_tempo(Voces, ns) #la funcion encuentra el tiempo maximo de termino de notas, es decir donde termina la pieza
    Tfinal = Array{Float64}(undef, ns)
    for i = 1:ns
        Tfinal[i] = maximum(Voces[i][:, 2])
    end
    return maximum(Tfinal)
end
"""
    get_onon_notes(s)
"""
function get_onon_notes(s)
    ini = s[s[:, 6].!=0, :] #getting the the starting places
    fin = s[s[:, 6].==0, :] #ending
    dn = unique(ini[:, 5]) #different notes.
    ndn = length(dn) #number of different notes
    posin_dn = Array{Vector}(undef, ndn) #initialize arrays for getting information of each note
    posfn_dn = Array{Vector}(undef, ndn)
    in_fi = Array{Matrix}(undef, ndn)
    for i = 1:ndn       #for each different note it gets when are they played
        tmp = ini[ini[:, 5].==dn[i], :][:, 2]
        tmp2 = unique(fin[fin[:, 5].==dn[i], :][:, 2])#added the unique to be sure they only finish "once", weird stuff was happening because of this.
        l1 = length(tmp)
        l2 = length(tmp2)
        if l1 != l2
            nt = min(l1, l2)
            note = [dn[i] for j = 1:nt]
            in_fi[i] = [tmp[1:nt] tmp2[1:nt] note]
        else
            note = [dn[i] for j = 1:length(tmp)]
            in_fi[i] = [tmp tmp2 note]
        end
    end
    nmat = vcat(in_fi...) #getting all notes in the same array
    return nmat[sortperm(nmat[:, 1]), :]  #returns the array sorted by appereance.
end

function isminor(s::Any)
    return occursin(r"^[a-z,#,/,\s]+$", s)
end

function isminor(k_seq::Array{Any,1})
    gkey = get_rank_freq(k_seq)[1, 1]
    return occursin(r"^[a-z,#,/,\s]+$", gkey)
end

"""
    funhar_seq(kseq, fun_key)
    Returns a sequence of Roman numerals given a reference for the Tonic Key.
"""
function funhar_seq(kseq, fun_key)
    fh_kseq = []
    if isminor(fun_key)
        nk = findfirst(x -> x == fun_key, nminor_keys)
    else
        nk = findfirst(x -> x == fun_key, nMajor_keys)
    end
    minork = circshift(nminor_keys, 13 - nk)
    majork = circshift(nMajor_keys, 13 - nk)
    for k = 1:length(kseq)
        if isminor(kseq[k])
            pos = findfirst(x -> x == kseq[k], minork) - 1
            push!(fh_kseq, Minor_RN[pos])
        else
            pos = findfirst(x -> x == kseq[k], majork) - 1
            push!(fh_kseq, Major_RN[pos])
        end
    end
    return fh_kseq
end
function circular_to_linear(S::Vector{Int})
    """
    Transform a set of notes into linear positions using circular statistics.

    Args:
        S: Vector of note values (natural numbers)
        
    Returns:
        L: Vector of linear positions (integers/floats)
    """

    # Step 1: Project to pitch classes once and count frequencies
    pitch_classes = [mod(n, 12) for n in S] 
    frequency = Dict{Int, Int}()
    for pc in pitch_classes
        frequency[pc] = get(frequency, pc, 0) + 1
    end

    # Step 2: Calculate circular weighted mean
    sum_cos = 0.0
    sum_sin = 0.0
    total_weight = length(S)  # More efficient than summing frequencies

    for (v, weight) in frequency
        angle = 2π * v / 12
        sum_cos += weight * cos(angle)
        sum_sin += weight * sin(angle)
    end

    X_bar = sum_cos / total_weight
    Y_bar = sum_sin / total_weight

    # Calculate circular mean
    circular_mean = mod(12 * atan(Y_bar, X_bar) / (2π), 12)

    # Step 3: Unwrap to linear positions (reuse computed pitch classes)
    L = Vector{Float64}(undef, length(S))  # Pre-allocate for better performance

    for (i, v) in enumerate(pitch_classes)
        distance = v - circular_mean

        if distance > 6
            L[i] = v - 12
        elseif distance < -6
            L[i] = v + 12
        else
            L[i] = v
        end
    end
    # Shift all values up by octave if any negative values exist
    if any(L .< 0)
        L = L .+ 12
    end

    return L
end
function cluster_notes(pitches)
    # Convert to pitch classes (0-11)
    pitch_classes = map(x -> mod(x, 12), pitches)

    # Initial conversion
    notes = get_cfpitch(pitch_classes)

    # Initial clustering by handling high/low split
    low_indices = findall(x -> x <= 6, notes)
    high_indices = findall(x -> x > 6, notes)

    if !isempty(high_indices) && !isempty(low_indices)
        # Shift the smaller group
        indices_to_shift = length(high_indices) <= length(low_indices) ? high_indices : low_indices
        for idx in indices_to_shift
            notes = shift_note_by_octave(notes, idx)
        end
    end

    # Iterative refinement
    best_notes = copy(notes)
    best_dispersion = calculate_dispersion(notes)

    # Set a maximum iteration count and convergence threshold
    max_iterations = 30
    convergence_threshold = 0.0001

    for _ in 1:max_iterations
        # Find the outlier more efficiently
        outlier_idx = find_outlier(notes)

        # Shift the outlier
        new_notes = shift_note_by_octave(notes, outlier_idx)

        # Calculate new dispersion
        new_dispersion = calculate_dispersion(new_notes)

        # Check if we improved
        if new_dispersion < best_dispersion - convergence_threshold
            best_notes = copy(new_notes)
            best_dispersion = new_dispersion
            notes = new_notes
        else
            # We've converged
            break
        end
    end

    return best_notes
end

function calculate_dispersion(notes)
    # Use variance as a measure of dispersion
    return Statistics.var(notes)
end

function find_outlier(notes)
    # Find the note that, when removed, minimizes the variance
    total_mean = Statistics.mean(notes)
    n = length(notes)

    # If we precompute sum, we can efficiently calculate mean without each element
    total_sum = sum(notes)

    max_deviation = -Inf
    outlier_idx = 1

    for i in 1:n
        # Mean without element i
        mean_without_i = (total_sum - notes[i]) / (n - 1)
        # How much this element deviates from mean of others
        deviation = abs(notes[i] - mean_without_i)

        if deviation > max_deviation
            max_deviation = deviation
            outlier_idx = i
        end
    end

    return outlier_idx
end

function shift_note_by_octave(notes, idx)
    new_notes = copy(notes)
    median_value = Statistics.median(notes)

    # Shift up or down by an octave
    if notes[idx] > median_value
        new_notes[idx] -= 12
    else
        new_notes[idx] += 12
    end

    return new_notes
end

function get_distance_ces(ce1, ce2)
    z_dif = ce2[3] - ce1[3]
    while abs(z_dif) > h_octav / 2 #translating over z to be in the same octave (same CE region)
        if z_dif > 0
            ce2[3] = ce2[3] - h_octav
        else
            ce2[3] = ce2[3] + h_octav
        end
        z_dif = ce2[3] - ce1[3]
    end
    return round(euclidean(ce1, ce2), digits=4)
end

"""
    get_local_lin_w(ptcs, factor)

    Returns a list of pitches and their respective weights from a sequence of pitches (ptcs)
    and a weighting factor (factor), the weighting gives more importance to lower notes
    within the 1-factor range.
"""
function get_local_lin_w(ptcs, factor)
    notas = sort(unique(ptcs))
    if length(notas) == 1
        n_we = []
        push!(n_we, 1)
    else
        p_min = notas[1]
        p_max = notas[end]
        tam = length(notas)
        n_we = zeros(tam)
        for i = 1:tam
            n_we[i] = 1 + (factor - 1) * (notas[i] - p_max) / (p_min - p_max)
        end
    end
    return notas, n_we
end

"""
    get_rank_freq(series)

    Returns a 2-Dimensional Array of symbols and their respective frequencies
    from a series of data. The output is ordered by rank (most to least frequent).
"""
function get_rank_freq(series)
    tam = length(series)
    M = Dict{Any,Int64}()
    for i = 1:tam
        M[series[i]] = get(M, series[i], 0) + 1
    end
    dist = sort(collect(M), by=tuple -> last(tuple), rev=true)
    rf = Array{Any}(undef, length(dist), 2)
    for i = 1:length(dist)
        rf[i, 1] = dist[i][1]
        rf[i, 2] = dist[i][2]
    end
    return rf
end