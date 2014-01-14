import re

class SearchTerm(object):
    pass

class Literal(SearchTerm):
    def __init__(self, term):
        self.term = term

    def search(self, string):
        return self.term in string

    def __str__(self):
        return self.term


class Regex(SearchTerm):
    def __init__(self, term):
        self.term = term
        self.match = re.compile(term)

    def search(self, string):
        return self.match.search(string) is not None

    def __str__(self):
        return self.term


def parse_search_term(string):
    search_term = string.strip()
    words = search_term.split()
    if len(words) >= 2 and words[0] == 'regex':
        return Regex(' '.join(words[1:]))
    elif len(search_term) > 0:
        return Literal(search_term)
    else:
        return None
