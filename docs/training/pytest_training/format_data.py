def format_data_for_display(people):
    """
    Format a list of people for display.

    Args:
        people (list): A list of dictionaries containing information about people.

    Returns:
        str: A formatted string representation of the people.
    """
    formatted_people = []
    for person in people:
        given_name = person.get("given_name", "Unknown")
        family_name = person.get("family_name", "Unknown")
        title = person.get("title", "Unknown")
        formatted_people.append(f"{given_name} {family_name}: {title}")
    return formatted_people