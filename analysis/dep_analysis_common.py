"""Functions common to dependence data analysis.
"""

# Interval around position that we search.
INTERVAL_SIZE = 300


def sort_seq_record_features_by_start(seq_record):
    seq_record.features = sorted(seq_record.features, key=lambda f: f.location.start)


def _is_position_in_interval(position, interval):
    return interval[0] <= position <= interval[1]


def _does_interval_overlap_feature(interval, feature):
    """Checks whether the given interval overlaps the feature's location.

    Args:
        interval: A two-tuple of integers (start, end).
        feature: A SeqFeature.

    Returns:
        A boolean indicating whether the features overlap.
    """
    interval_start = interval[0]
    interval_end = interval[1]

    if feature.location.start == interval_start:
        # >>>>>>>>>
        # (....)
        return interval_end - interval_start > 0

    elif feature.location.start < interval_start:
        # >>>>>>>>>
        #    (..........)
        return feature.location.end > interval_start

    else:
        #      >>>>>>>>>
        # (........)
        return feature.location.start < interval_end


def find_features_overlapped_by_position(position, features_sorted_seq_record):
    features_overlapped = []
    for feature in features_sorted_seq_record.features:
        if feature.location.start > position:
            # No need to look further since features are sorted by
            # start position.
            break
        feature_interval = (feature.location.start, feature.location.end)
        if _is_position_in_interval(position, feature_interval):
            features_overlapped.append(feature)
    return features_overlapped


def find_features_overlapped_by_interval(interval, features_sorted_seq_record):
    features_overlapped = []
    for feature in features_sorted_seq_record.features:
        if feature.location.start > interval[1]:
            # No need to look further since features are sorted by
            # start position.
            break
        if _does_interval_overlap_feature(interval, feature):
            features_overlapped.append(feature)
    return features_overlapped


IGNORE_FEATURE_TYPES = [
    'source',
    'repeat_region',
    'mobile_element',
    'variation',
    'STS'
]

def filter_features(features):
    """Remove duplicates and other unwanted features.
    """
    # Remove ignored types
    filtered = [f for f in features if not f.type in IGNORE_FEATURE_TYPES] 
    
    return filtered


def bucket_nearby_features(features, position):
    """Sorts into upstream and downstream, filtering out those that are overlapped.
    """
    upstream_of = []
    downstream_of = []
    for feature in features:
        # xxxxxxxx .
        if feature.location.end < position:
            # >>>>>>>> .
            if feature.strand == 1:
                downstream_of.append(feature)
            # <<<<<<<< .
            else:
                upstream_of.append(feature)
        # . xxxxxxxxx
        else:
            assert feature.location.start > position
            # . >>>>>>>>>
            if feature.strand == 1:
                upstream_of.append(feature)
            # . <<<<<<<<<
            else:
                downstream_of.append(feature)
        
    return {
        'snv_is_upstream_of': upstream_of,
        'snv_is_downstream_of': downstream_of
    }


def get_feature_label(feature):
    try_name_qualifier_keys = ['label', 'gene', 'note', 'operon', 'product']
    label = None
    for key in try_name_qualifier_keys:
        if key in feature.qualifiers:
            label = feature.qualifiers[key][0]
            break
    if not label:
        print feature, feature.qualifiers, key
        assert False, "Label not found"

    # Truncate label if necessary.
    if len(label) > 20:
        label = label[:10] + '...'

    # Add type in parens.
    label = label + ' (' + feature.type + ')'

    return label
