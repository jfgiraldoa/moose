#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

"""Defines Renderer objects that convert AST (from Reader) into an output format."""
import os
import re
import logging
import traceback
import codecs
import shutil
import anytree

import MooseDocs
from MooseDocs import common
from MooseDocs.common import exceptions, mixins
from MooseDocs.tree import html, latex, pages, tokens

LOG = logging.getLogger(__name__)

def create_directory(location):
    """Helper for creating a directory."""
    with MooseDocs.base.Translator.LOCK:
        dirname = os.path.dirname(location)
        if dirname and not os.path.isdir(dirname):
            LOG.debug('CREATE DIR %s', dirname)
            os.makedirs(dirname)

class Renderer(mixins.ConfigObject, mixins.ComponentObject):
    """
    Base renderer for converting AST to an output format.
    """

    #:[str] The name of the method to call on RendererComponent objects.
    METHOD = None

    def __init__(self, **kwargs):
        mixins.ConfigObject.__init__(self, **kwargs)
        mixins.ComponentObject.__init__(self)
        self.__functions = dict()  # functions on the RenderComponent to call

    def add(self, name, component):
        """
        Associate a RenderComponent object with a token type.

        Inputs:
            name[str]: The token name (e.g., "String") to associate with the supplied component.
            compoment[RenderComponent]: The component to execute with the associated token type.
        """
        if MooseDocs.LOG_LEVEL == logging.DEBUG:
            common.check_type("name", name, unicode)
            common.check_type("component", component, MooseDocs.base.components.RenderComponent)
        component.setRenderer(self)
        self.addComponent(component)
        self.__functions[name] = self._method(component)

    def getRoot(self):
        """
        Return the rendered content root node.

        Called by the Translator prior to beginning rendering.
        """
        raise NotImplementedError()

    def render(self, parent, token, page):
        """
        Convert the AST defined in the token input into a output node of parent.

        Inputs:
            ast[tree.token]: The AST to convert.
            parent[tree.base.NodeBase]: A tree object that the AST shall be converted to.
        """
        try:
            func = self.__getFunction(token)
            el = func(parent, token, page) if func else parent

        except Exception as e: #pylint: disable=broad-except
            el = None
            if token.info is not None:
                line = token.info.line
                src = token.info[0]
            else:
                line = None
                src = ''
            msg = common.report_error(e.message,
                                      page.source,
                                      line,
                                      src,
                                      traceback.format_exc(),
                                      u'RENDER ERROR')
            with MooseDocs.base.translators.Translator.LOCK:
                LOG.error(msg)

        if el is not None:
            for child in token.children:
                self.render(el, child, page)

    def preRender(self, result, page, meta):
        """
        Called by Translator prior to rendereing.

        Inputs:
            result[tree.base.NodeBase]: The root node of the result tree.
        """
        pass

    def postRender(self, result, page, meta):
        """
        Called by Translator after rendereing.

        Inputs:
            result[tree.base.NodeBase]: The root node of the result tree.
        """
        pass

    def write(self, page, result=None):
        """
        Write the supplied results using to the destination defined by the page.

        This is called by the Tranlator object.
        """
        if isinstance(page, pages.Source):
            create_directory(page.destination)
            LOG.debug('WRITE %s-->%s', page.source, page.destination)
            with codecs.open(page.destination, 'w', encoding='utf-8') as fid:
                fid.write(result.write())

        elif isinstance(page, pages.File):
            create_directory(page.destination)
            LOG.debug('COPY: %s-->%s', page.source, page.destination)
            if not os.path.exists(page.source):
                LOG.error('Unknown file: %s', page.source)
            else:
                shutil.copyfile(page.source, page.destination)

        elif isinstance(page, pages.Directory):
            create_directory(page.destination)

        else:
            LOG.error('Unknown Node type: %s', type(page))

    def _method(self, component):
        """
        Return the desired method to call on the RenderComponent object.

        Inputs:
            component[RenderComponent]: Object to use for locating desired method for renderering.
        """
        if self.METHOD is None:
            msg = "The Reader class of type {} must define the METHOD class member."
            raise exceptions.MooseDocsException(msg, type(self))
        elif not hasattr(component, self.METHOD):
            msg = "The component object {} does not have a {} method."
            raise exceptions.MooseDocsException(msg, type(component), self.METHOD)
        return getattr(component, self.METHOD)

    def __getFunction(self, token):
        """
        Return the desired function for the supplied token object.

        Inputs:
            token[tree.token]: token for which the associated RenderComponent function is desired.
        """
        return self.__functions.get(token.name, None)

class HTMLRenderer(Renderer):
    """
    Converts AST into HTML.
    """
    METHOD = 'createHTML'
    EXTENSION = '.html'

    @staticmethod
    def defaultConfig():
        """
        Return the default configuration.
        """
        config = Renderer.defaultConfig()
        config['google_analytics'] = (False, "Enable Google Analytics.")
        config['favicon'] = (None, "The location of the website favicon.")
        return config

    def __init__(self, *args, **kwargs):
        Renderer.__init__(self, *args, **kwargs)
        self.__head_javascript = common.Storage()
        self.__javascript = common.Storage()

        self.__css = common.Storage()

        if self.get('google_analytics', False):
            self.addJavaScript('google_analytics', 'js/google_analytics.js')

    def getRoot(self):
        """Return the result node for inserting rendered html nodes."""
        root = html.Tag(None, '!DOCTYPE html', close=False)
        html.Tag(root, 'head')
        return html.Tag(root, 'body')

    def addJavaScript(self, key, filename, location='_end', head=False, **kwargs):
        """Add a javascript dependency."""
        if head:
            self.__head_javascript.add(key, (filename, kwargs), location)
        else:
            self.__javascript.add(key, (filename, kwargs), location)

    def addCSS(self, key, filename, location='_end', **kwargs):
        """Add a CSS dependency."""
        self.__css.add(key, (filename, kwargs), location)

    def postRender(self, result, page, meta):
        """Insert CSS/JS dependencies into html node tree."""
        root = result.root

        def rel(path):
            """Helper to create relative paths for js/css dependencies."""
            return os.path.relpath(path, os.path.dirname(page.local))

        head = anytree.search.find_by_attr(root, 'head')
        body = anytree.search.find_by_attr(root, 'body')

        favicon = self.get('favicon')
        if favicon:
            html.Tag(head, 'link', rel="icon", type="image/x-icon", href=rel(favicon), \
                     sizes="16x16 32x32 64x64 128x128")

        for name, kwargs in self.__css:
            html.Tag(head, 'link', href=rel(name), type="text/css", rel="stylesheet", **kwargs)

        for name, kwargs in self.__head_javascript:
            html.Tag(head, 'script', type="text/javascript", src=rel(name), **kwargs)

        for name, kwargs in self.__javascript:
            html.Tag(body.parent, 'script', type="text/javascript", src=rel(name), **kwargs)

class MaterializeRenderer(HTMLRenderer):
    """
    Convert AST into HTML using the materialize javascript library (http://materializecss.com).
    """
    METHOD = 'createMaterialize'

    @staticmethod
    def defaultConfig():
        """
        Return the default configuration.
        """
        config = HTMLRenderer.defaultConfig()
        config['breadcrumbs'] = (True, "Toggle for the breadcrumb links at the top of page.")
        config['sections'] = (True, "Group heading content into <section> tags.")
        config['collapsible-sections'] = ([None, None, None, None, None, None],
                                          "Collapsible setting for the six heading level " \
                                          "sections, possible values include None, 'open', and " \
                                          "'close'. Each indicates if the associated section " \
                                          "should be collapsible, if so should it be open or " \
                                          "closed initially. The 'sections' setting must be " \
                                          "True for this to operate.")
        config['navigation'] = (None, "Top bar website navigation items.")
        config['repo'] = (None, "The source code repository.")
        config['name'] = (None, "The name of the website (e.g., MOOSE)")
        config['home'] = ('/', "The homepage for the website.")
        config['scrollspy'] = (True, "Enable/disable the scrolling table of contents.")
        config['search'] = (True, "Enable/disable the search bar.")
        return config

    def __init__(self, *args, **kwargs):
        HTMLRenderer.__init__(self, *args, **kwargs)
        self.__navigation = None # cache for navigation pages
        self.__index = False     # page index created

        self.addCSS('materialize', "contrib/materialize/materialize.min.css",
                    media="screen,projection")
        self.addCSS('prism', "contrib/prism/prism.min.css")
        self.addCSS('moose', "css/moose.css")

        self.addJavaScript('jquery', "contrib/jquery/jquery.min.js")
        self.addJavaScript('materialize', "contrib/materialize/materialize.min.js")
        self.addJavaScript('clipboard', "contrib/clipboard/clipboard.min.js")
        self.addJavaScript('prism', "contrib/prism/prism.min.js")
        self.addJavaScript('init', "js/init.js")

    def update(self, **kwargs):
        """
        Update the default configuration with the supplied values. This is an override of the
        ConfigObject method and is simply modified here to the check the type of a configuration
        item.
        """
        HTMLRenderer.update(self, **kwargs)

    def getRoot(self):
        body = HTMLRenderer.getRoot(self)

        wrap = html.Tag(body, 'div', class_='page-wrap')
        html.Tag(wrap, 'header')

        main = html.Tag(wrap, 'main', class_='main')
        container = html.Tag(main, 'div', class_="container")

        row = html.Tag(container, 'div', class_="row")
        col = html.Tag(row, 'div', class_="moose-content")

        return col

    def _method(self, component):
        """
        Fallback to the HTMLRenderer method if the MaterializeRenderer method is not located.

        Inputs:
            component[RenderComponent]: Object to use for locating desired method for renderering.
        """
        if hasattr(component, self.METHOD):
            return getattr(component, self.METHOD)
        elif hasattr(component, HTMLRenderer.METHOD):
            return getattr(component, HTMLRenderer.METHOD)
        else:
            msg = "The component object {} does not have a {} method."
            raise exceptions.MooseDocsException(msg, type(component), self.METHOD)

class LatexRenderer(Renderer):
    """
    Renderer for converting AST to LaTeX.
    """
    METHOD = 'createLatex'
    EXTENSION = '.tex'

    def __init__(self, *args, **kwargs):
        self._packages = dict()
        self._preamble = list()
        self._commands = dict()
        Renderer.__init__(self, *args, **kwargs)

    def getRoot(self):
        """
        Return LaTeX root node.
        """
        return latex.LatexBase(None, None)

    def addNewCommand(self, cmd, content):
        """
        Add a NewDocumentCommand to latex preamble.
        """
        num = 0
        for match in re.finditer(r'#(?P<num>[0-9]+)', content):
            num = max(num, int(match.group('num')))

        args = [latex.Brace(string=cmd, escape=False), latex.Brace(string='m'*num)]
        self._commands[cmd] = latex.Command(None, 'NewDocumentCommand', args=args, escape=False,
                                            string=content, start='\n')

    def getNewCommands(self):
        """Return the dict of new commands."""
        return self._commands

    def addPackage(self, pkg, *args, **kwargs):
        """
        Add a LaTeX package to the list of packages for rendering (see pdf.py)
        """
        self._packages[pkg] = (args, kwargs)

    def getPackages(self):
        """Return the set of packages and settings."""
        return self._packages

    def addPreamble(self, node):
        """
        Add a string to the preamble (see pdf.py).
        """
        self._preamble.append(node)

    def getPreamble(self):
        """Return the list of preamble strings."""
        return self._preamble

class RevealRenderer(HTMLRenderer):
    """
    Convert AST into HTML using the materialize javascript library (http://materializecss.com).
    """
    METHOD = 'createReveal'

    @staticmethod
    def defaultConfig():
        """
        Return the default configuration.
        """
        config = HTMLRenderer.defaultConfig()
        config['theme'] = ('simple', "The CSS theme to use (simple).")
        return config

    def __init__(self, *args, **kwargs):
        HTMLRenderer.__init__(self, *args, **kwargs)
        self.addCSS('reveal', "contrib/reveal/reveal.css")

        self.addCSS('reveal_theme', "contrib/reveal/{}.css".format(self.get('theme')), id_="theme")
        self.addCSS('prism', "contrib/prism/prism.min.css")
        self.addCSS('reveal_css', "css/reveal_moose.css", location='>reveal_theme')

        self.addJavaScript('reveal', "contrib/reveal/reveal.js")
        self.addJavaScript('prism', "contrib/prism/prism.min.js")
        self.addJavaScript('notes', "contrib/reveal/notes.js")
        self.addJavaScript('reveal_init', "js/reveal_init.js")

    def getRoot(self):
        body = HTMLRenderer.getRoot(self)
        div = html.Tag(body, 'div', class_='reveal')
        slides = html.Tag(div, 'div', class_='slides')
        return slides#html.Tag(slides, 'section')

    def _method(self, component):
        if hasattr(component, self.METHOD):
            return getattr(component, self.METHOD)
        elif hasattr(component, HTMLRenderer.METHOD):
            return getattr(component, HTMLRenderer.METHOD)
        else:
            msg = "The component object {} does not have a {} method."
            raise exceptions.MooseDocsException(msg, type(component), self.METHOD)


class JSONRenderer(Renderer):
    """
    Render the AST as a JSON file.
    """
    METHOD = 'dummy' # The AST can write JSON directly.
    EXTENSION = '.json'

    def getRoot(self): #pylint: disable=unused-argument
        """
        Return LaTeX root node.
        """
        return tokens.Token()

    def _method(self, component):
        """
        Do not call any methods on the component, just create a new tree for writing JSON.
        """
        return self.__function

    @staticmethod
    def __function(parent, token, page):
        """Replacement for Component function (see _method)."""
        token.parent = parent
        return token
