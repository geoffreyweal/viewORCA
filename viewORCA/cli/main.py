#from viewORCA import __version__
#from viewORCA.check_prerequisites import check_prerequisites
import argparse, textwrap
from importlib import import_module

class CLIError(Exception):
    """Error for CLI commands.

    A subcommand may raise this.  The message will be forwarded to
    the error() method of the argument parser."""


# Important: Following any change to command-line parameters, use
# python3 -m ase.cli.completion to update autocompletion.
commands = [
    ('multi',    'viewORCA.Programs.get_multiplicity'),
    ('opt',      'viewORCA.Programs.viewOPT'),
    ('scan',     'viewORCA.Programs.viewSCAN'),
    ('neb',      'viewORCA.Programs.viewNEB'),
    ('neb_snap', 'viewORCA.Programs.viewNEB_snapshot'),
    ('irc',      'viewORCA.Programs.viewIRC'),
]

def main(prog='viewORCA', description='Program for viewing jobs from ORCA', commands=commands, hook=None, args=None):
    """
    This method is designed to allow the user to view jobs performed by ORCA to study organic mechanisms. 
    """
    #check_prerequisites(__version__)

    parser = argparse.ArgumentParser(prog=prog, description=description, formatter_class=Formatter)
    parser.add_argument('-T', '--traceback', action='store_true')
    subparsers = parser.add_subparsers(title='Sub-commands', dest='command')
    subparser = subparsers.add_parser('help', description='Help', help='Help for sub-command.')
    subparser.add_argument('helpcommand', nargs='?', metavar='sub-command', help='Provide help for sub-command.')

    functions = {}
    parsers = {}
    for command, module_name in commands:
        cmd = import_module(module_name).CLICommand
        docstring = cmd.__doc__
        if docstring is None:
            # Backwards compatibility with GPAW
            short = cmd.short_description
            long = getattr(cmd, 'description', short)
        else:
            parts = docstring.split('\n', 1)
            if len(parts) == 1:
                short = docstring
                long = docstring
            else:
                short, body = parts
                long = short + '\n' + textwrap.dedent(body)
        subparser = subparsers.add_parser(
            command,
            formatter_class=Formatter,
            help=short,
            description=long)
        cmd.add_arguments(subparser)
        functions[command] = cmd.run
        parsers[command] = subparser

    if hook:
        args = hook(parser, args)
    else:
        args = parser.parse_args(args)

    if args.command == 'help':
        if args.helpcommand is None:
            parser.print_help()
        else:
            parsers[args.helpcommand].print_help()
    elif args.command is None:
        parser.print_usage()
    else:
        f = functions[args.command]
        try:
            if f.__code__.co_argcount == 1:
                f(args)
            else:
                f(args, parsers[args.command])
        except KeyboardInterrupt:
            pass
        except CLIError as x:
            parser.error(x)
        except Exception as x:
            if args.traceback:
                raise
            else:
                l1 = f'{x.__class__.__name__}: {x}\n'
                l2 = ('To get a full traceback, use: {} -T {} ...'
                      .format(prog, args.command))
                parser.error(l1 + l2)


class Formatter(argparse.HelpFormatter):
    """Improved help formatter."""

    def _fill_text(self, text, width, indent):
        assert indent == ''
        out = ''
        blocks = text.split('\n\n')
        for block in blocks:
            if block[0] == '*':
                # List items:
                for item in block[2:].split('\n* '):
                    out += textwrap.fill(item,
                                         width=width - 2,
                                         initial_indent='* ',
                                         subsequent_indent='  ') + '\n'
            elif block[0] == ' ':
                # Indented literal block:
                out += block + '\n'
            else:
                # Block of text:
                out += textwrap.fill(block, width=width) + '\n'
            out += '\n'
        return out[:-1]
