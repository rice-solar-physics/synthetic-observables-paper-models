"""
Parse a LaTeX file, make figure paths relative, and remove other PythonTeX stuff
"""
import os
import parse
import click


@click.group()
def cli():
    click.echo('Fixing LaTeX document...')


@cli.command()
@click.option('--input-file', '-i', multiple=True, help='TeX file to parse')
@click.option('--output-file', '-o', multiple=True, help='Path to fixed TeX file')
def fix_figure_paths(input_file, output_file):
    with open(input_file[0], 'r') as f:
        lines = f.readlines()

    for i, l in enumerate(lines):
        if 'includegraphics' in l:
            r = parse.parse('\includegraphics[{}]{{{}}}', l.strip())
            _, fig_filename = os.path.split(r[1])
            new_line = f'   \includegraphics[{r[0]}]{{{fig_filename}}}\n'
            lines[i] = new_line

    with open(output_file[0], 'w') as f:
        f.write(''.join(lines))


@cli.command()
@click.option('--input-file', '-i', multiple=True, help='TeX file to parse')
@click.option('--output-file', '-o', multiple=True, help='Path to fixed TeX file')
def remove_pythontex_block(input_file, output_file):
    with open(input_file[0], 'r') as f:
        lines = f.readlines()

    for i, l in enumerate(lines):
        if 'begin pythontex bug fix' in l.lower():
            i_start = i
        if 'end pythontex bug fix' in l.lower():
            i_end = i

    with open(output_file[0], 'w') as f:
        f.write(''.join(lines[:i_start] + lines[i_end+1:]))


if __name__ == '__main__':
    cli()
