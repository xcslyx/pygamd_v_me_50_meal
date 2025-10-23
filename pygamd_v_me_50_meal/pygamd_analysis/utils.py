def str2value(value):
    if value.lower() in ('yes', 'y', 'true', 'ture', 't', '1'):
        return True
    elif value.lower() in ('no', 'n', 'false', 'f', '0'):
        return False
    elif value == "avg":
        return "avg"
    elif value == "unset":
        return None
    else:
        return value
