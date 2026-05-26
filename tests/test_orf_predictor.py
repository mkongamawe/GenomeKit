from genomekit.modules.orf_predictor import find_orfs


def test_find_orfs_single_orf():
    frames = [
        ["ATG", "CCC", "TAG"],  # Frame 1
        ["CCC", "ATG", "AAA"],  # Frame 2
        ["GGG", "TTT", "CCC"],  # Frame 3
    ]

    result = find_orfs(frames)

    expected = [(1, "ATGCCCTAG")]

    assert result == expected


def test_find_orfs_multiple_orfs():
    frames = [
        ["ATG", "AAA", "TAA", "ATG", "CCC", "TAG"],
        ["ATG", "GGG", "TGA"],
    ]

    result = find_orfs(frames)

    expected = [
        (1, "ATGAAATAA"),
        (1, "ATGCCCTAG"),
        (2, "ATGGGGTGA"),
    ]

    assert result == expected


def test_find_orfs_no_stop_codon():
    frames = [["ATG", "CCC", "AAA"]]

    result = find_orfs(frames)

    assert result == []


def test_find_orfs_no_start_codon():
    frames = [["CCC", "AAA", "TAG"]]

    result = find_orfs(frames)

    assert result == []


def test_find_orfs_empty_frames():
    frames = []

    result = find_orfs(frames)

    assert result == []
