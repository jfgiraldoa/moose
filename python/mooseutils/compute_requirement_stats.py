import os
from mooseutils import colorText, git_ls_files, git_root_dir
try:
    from hit_load import hit_load
except:
    raise ImportError('hit package failed to load, but it is required.')

class SQAStats(object):
    """Data wrapper for SQA statistic information."""
    def __init__(self, title):
        self.title = title
        self.files = 0
        self.files_with_requirement = 0
        self.tests = 0
        self.tests_with_requirement = 0

    @property
    def complete(self):
        return float(self.tests_with_requirement)/float(self.tests) if self.tests else 0

    def __iadd__(self, other):
        self.files += other.files
        self.files_with_requirement += other.files_with_requirement
        self.tests += other.tests
        self.tests_with_requirement += other.tests_with_requirement
        return self

    def __str__(self):
        """Display requirement stats."""
        out = []
        out.append(colorText(self.title, 'LIGHT_YELLOW'))
        out.append('                   Complete: {:2.1f}%'.format(self.complete*100))
        out.append('      Total Number of Files: {}'.format(self.files))
        out.append('    Files with Requirements: {}'.format(self.files_with_requirement))
        out.append('      Total Number of Tests: {}'.format(self.tests))
        out.append('    Tests with Requirements: {}'.format(self.tests_with_requirement))
        return '\n'.join(out)

def compute_requirement_stats(location, specs=['tests'], working_dir=None, show=True):
    """
    Report requirement statistics for the test spec files with the supplied location.

    Inputs:
        location: Path to directory contain test specifications, the supplied path should be
                  relative to the cwd input.
        specs: The filename(s) to consider
        working_dir: The working directory, if not supplied the root directory of the repository is used

    """
    working_dir = git_root_dir(os.getcwd()) if working_dir is None else working_dir
    data = SQAStats(location)
    location = os.path.join(working_dir, location)
    for filename in git_ls_files(location):
        if not os.path.isfile(filename):
            continue
        fname = os.path.basename(filename)
        if fname in specs:
            root = hit_load(filename)
            has_requirement = False
            for child in root.children[0]:
                data.tests += 1
                if child.get('requirement', None):
                    has_requirement = True
                    data.tests_with_requirement += 1

            data.files += 1
            if has_requirement:
                data.files_with_requirement += 1

    if show:
        print(data)
    return data
