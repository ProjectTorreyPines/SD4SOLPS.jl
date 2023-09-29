using SD4SOLPS: SD4SOLPS
using SOLPS2IMAS: SOLPS2IMAS

eqdsk_file = "g002296.00200"
sample_paths =
    [
        splitdir(pathof(SD4SOLPS))[1] * "/../sample",
        splitdir(pathof(SOLPS2IMAS))[1] * "/../samples",
    ]
test_dir = mktempdir()
ids = SD4SOLPS.preparation(eqdsk_file, sample_paths...; filename=test_dir * "/output")
# core_profile_2d(dd, prof_time_idx, eq_time_idx, quantity, r, z)
SD4SOLPS.core_profile_2d(ids, 1, 1, "electrons.density", 3.4, 7.0)