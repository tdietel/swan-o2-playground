from IPython.core.magic import register_cell_magic, no_var_expand
# import string
# import re
import argparse
import pathlib
import jinja2


@register_cell_magic
@no_var_expand
def genfile(line, cell):
    """ Generate file 
    
    This magic is an extension of the %%writefile magic. The additional functionality 
    includes:
      - variable expansion wih Jinja2
      - set file permisssions with -m | --mode flag"""

    # PARSE COMMAND LINE ARGUMENTS
    parser = argparse.ArgumentParser(description='Generate file magic')
    parser.add_argument('outfile', type=str, help='output file')
    parser.add_argument('-m', '--mode', help='file permissions')
    args = parser.parse_args(line.split())

    # FORMAT THE CONTENT
    
    # We ask get_python() to give us all globals with ev("globals()"), which 
    # we then provide to the Jinja2 template
    content = jinja2.Template(cell).render(**get_ipython().ev("globals()"))
    
    # Replace leading spaces with tabs
    # TODO: it would be better to make this configurable with
    # an argparse option
    content = content.replace("\n    ", "\n\t")

    # TODO: this should check if the file content changed to
    # prevent updating the timestamp if nothing changed.
    # This will help with Makefile dependencies.

    # GENERATE THE FILE    
    filename = pathlib.Path(args.outfile)
    with open(filename, "w") as out:
        out.write(content)
    if args.mode is not None:
        filename.chmod(int(args.mode, 8))

    print(f"Generated file '{filename}'")

    return
