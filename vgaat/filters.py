def filter_by_alpha(result_groups: dict[str, list[dict]], alpha=0.05):
    filtered_result_groups = {}
    for identifier, result_group in result_groups.items():
        test_alias, viral_strand, regions, synonymity, fishers = identifier
        filtered_result_data = []
        title, desc, result_data = result_group
        for info in result_data:
            if (not isinstance(info["p"], str)) and info["p"] < float(alpha):
                filtered_result_data.append(dict(info))
        if len(filtered_result_data) > 0:
            filtered_result_groups[identifier] = (title, desc, filtered_result_data)
    return filtered_result_groups
