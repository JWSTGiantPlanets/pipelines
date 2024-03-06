import argparse


class BooleanOrBothAction(argparse.Action):
    """
    Based on BooleanOptionalAction, but with the ability to set the value to 'both'.
    """

    def __init__(
        self,
        option_strings,
        dest,
        default=None,
        type=None,
        choices=None,
        required=False,
        help=None,
        metavar=None,
    ):

        _option_strings = []
        for option_string in option_strings:
            _option_strings.append(option_string)

            if option_string.startswith('--'):
                option_string_alt = '--no-' + option_string[2:]
                _option_strings.append(option_string_alt)

                option_string_alt = '--both-' + option_string[2:]
                _option_strings.append(option_string_alt)

        super().__init__(
            option_strings=_option_strings,
            dest=dest,
            nargs=0,
            default=default,
            type=type,
            choices=choices,
            required=required,
            help=help,
            metavar=metavar,
        )

    def __call__(
        self,
        parser,
        namespace,
        values,
        option_string=None,
    ):
        if option_string in self.option_strings:
            if option_string.startswith('--both-'):
                value = 'both'
            else:
                value = not option_string.startswith('--no-')
            setattr(namespace, self.dest, value)

    def format_usage(self):
        return ' | '.join(self.option_strings)
