"""
Parse a LaTeX file and make figure paths relative
"""
import os
import parse
import click


@click.command()
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


if __name__ == '__main__':
    fix_figure_paths()
