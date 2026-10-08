---
search:
  exclude: true
---

# How an offline document navigates with no server

Research note for issue #397 (parent #395). Everything below was read or measured on
**2026-10-07** unless a line says otherwise. Browser measurements are Google Chrome 154 on
macOS 15 (`navigator.userAgent`: `Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)
AppleWebKit/537.36 (KHTML, like Gecko) Chrome/154.0.0.0 Safari/537.36`).

`docs/adr/0002-editable-protocols.md` and `CONTEXT.md` both commit a protocol to one
self-contained HTML page. #393's point 3 wants a project split across pages with a nav bar.
This note gathers what the three routes actually cost, what published protocol formats do, how
the layout #393 names behaves, and what printing means. **It decides nothing** — #399 decides.

## 1. What is on the page today

Measured on this repository, `docs/examples/`:

| | `pUC19-GFP` (Golden Gate) | `ap1-library` (IGGA) |
| --- | --- | --- |
| `protocol.html` | 63,790 bytes | 364,517 bytes |
| `protocol.json` | 38,177 bytes | 236,428 bytes |
| `<script src=>` to anywhere | 0 | 0 |
| `<link>` to anywhere | 0 | 0 |
| `id=` anchors in the page | — | 73 |
| Printed pages, US Letter | **16** | **66** |

Page counts are Chrome 154 headless, `--print-to-pdf --no-pdf-header-footer`, default margins,
counted with `pypdf`. The pages are one file each: nothing is fetched, so the file opens from a
mail attachment with no network.

Navigation inside that file is one flat list. `render.py:139-143` writes a `<nav class="toc">`
of one `<a href="#step-N">` per step — 42 of them for AP-1 — above the first step, and nothing
else. There are no sections and no second level.

The page already prints deliberately. `protocol.css:209-241` is one `@media print` block: it
forces a light palette, hides `.toolbar`, `.toc` and the copy buttons, and **opens every closed
disclosure** —

```css
.oligo-checks:not([open]) > ul { display: block; }
.oligo-checks::details-content { content-visibility: visible; }
```

with the comment *"Paper has no toggle, so the reasons print with the sheet."* That is route 3's
printing problem already met once, in this repository, and solved by hand per element.

## 2. The three routes, under `file://`

The measurement below is a probe page opened from disk in Chrome 154. Its four files sit in one
directory; it reports what each API does.

| Probe | Result under `file://` |
| --- | --- |
| `window.origin` | `"null"` |
| `window.isSecureContext` | `true` |
| `<script src="./c.js">` (classic) | ran |
| `<script type="module" src="./m.js">` | **did not run** |
| `fetch('./data.json')` | **threw** `TypeError: Failed to fetch` |
| `XMLHttpRequest` on `./data.json`, same directory | **threw** `Failed to load 'file:///…/data.json'` |
| `new Worker('./w.js')` | **threw** `Script at 'file:///…/w.js' cannot be accessed from origin 'null'` |
| `history.pushState({}, '', '#page2')` | ok — URL became `…/index.html#page2` |
| `history.pushState({}, '', 'second.html')` | **threw** `SecurityError: … cannot be created in a document with origin 'null'` |
| `localStorage.setItem` | ok |
| `<a href="./second.html">` clicked | navigated |
| Back button after that navigation | returned to `…/index.html#page2` |

The `SecurityError` text is worth keeping whole, because it is the finding:

> Failed to execute 'pushState' on 'History': A history state object with URL
> 'file:///…/second.html' cannot be created in a document with origin 'null' and URL
> 'file:///…/index.html#page2'.

### Why, in the specs

**The origin is the cause.** MDN: *"Modern browsers usually treat the origin of files loaded
using the `file:///` scheme as opaque origins. What this means is that if a file includes other
files from the same folder (say), they are not assumed to come from the same origin, and may
trigger CORS errors"*
(<https://developer.mozilla.org/en-US/docs/Web/Security/Same-origin_policy>). The CORS error page
is more specific: *"CORS requests may only use the HTTP or HTTPS URL scheme […] This often occurs
if the URL specifies a local file, using the `file:///` scheme"*, and *"Many browsers, including
Firefox and Chrome, now treat all local files as having opaque origins (by default)"*
(<https://developer.mozilla.org/en-US/docs/Web/HTTP/Guides/CORS/Errors/CORSRequestNotHttp>). That
page names the APIs CORS governs — `fetch()`, `XMLHttpRequest`, `@font-face`, WebGL textures and
XSL stylesheets — and classic `<script src>`, `<link rel=stylesheet>` and `<img>` are not among
them, which is what the probe measured.

Firefox dates the change: CVE-2019-11730, *"Same-origin policy treats all files in a directory as
having the same-origin"*, fixed in Firefox 68
(<https://www.mozilla.org/en-US/security/advisories/mfsa2019-21/>). Chromium's switch table says
the same in a comment on `kAllowFileAccessFromFiles`: *"By default, file:// URIs cannot read other
file:// URIs. This is an override for developers who need the old behavior for testing"*
(<https://chromium.googlesource.com/chromium/src/+/main/content/public/common/content_switches.cc>)
— a testing flag, not a deployment mode.

The Fetch standard simply declines the case. Its scheme fetch step reads *"For now, unfortunate
as it is, `file:` URLs are left as an exercise for the reader"*, and the fall-through is *"Return
a network error"* (<https://fetch.spec.whatwg.org/>).

**Module scripts are the sharpest edge.** MDN, twice on one page: *"if you try to load the HTML
file locally (i.e., with a `file://` URL), you'll run into CORS errors due to JavaScript module
security requirements. You need to do your testing through a server"*
(<https://developer.mozilla.org/en-US/docs/Web/JavaScript/Guide/Modules>). Measured above: the
module never ran; the classic script did.

**`pushState` is governed by the URL, not the origin, and the HTML standard carves out `file:`
by hand.** Its "can have its URL rewritten" algorithm returns false when the two URLs *"differ in
their scheme, username, password, host, or port components"*, returns true for an HTTP(S) scheme,
and for `file:` says: *"If targetURL and documentURL differ in their path component, then return
false […] Return true"*, with the note *"Differences in query and fragment are allowed for
'file:' URLs"* (<https://html.spec.whatwg.org/multipage/nav-history-apis.html>). So a router that
pushes `#step-3` or `?step=3` survives on disk; one that pushes `step-3.html` throws, exactly as
measured.

### What route 3 costs, in find-in-page

Hiding sections inside one file costs the browser's own search, and the rules are per mechanism.

- `content-visibility: hidden` — *"The skipped contents must not be accessible to user-agent
  features, such as find-in-page […] This is similar to giving the contents `display: none`."*
  `content-visibility: auto` — *"the skipped contents must still be available as normal to
  user-agent features such as find-in-page"*
  (<https://developer.mozilla.org/en-US/docs/Web/CSS/content-visibility>).
- `hidden="until-found"` is the escape hatch: *"the element is hidden but its content will be
  accessible to the browser's 'Find in page' feature or to fragment navigation"*, implemented
  with `content-visibility: hidden`. Its trap is in the same page: *"If the element in the hidden
  until found state has a `display` value of `none`, `contents`, or `inline`, then the element
  will not be revealed by 'Find in page' or fragment navigation."*
  (<https://developer.mozilla.org/en-US/docs/Web/HTML/Reference/Global_attributes/hidden>).
- The HTML standard makes this normative in its find-in-page section: elements *"with the `hidden`
  attribute in the Hidden Until Found state should have their skipped contents become
  accessible"*, and `details` elements *"which do not have their open attribute set should have
  the skipped contents of their second slot"*, restored synchronously so it is *"not observable to
  users or to author code"* (<https://html.spec.whatwg.org/multipage/interaction.html>).
- Chrome shipped the `<details>` half in 97 — *"Closed `<details>` elements are now searchable and
  can be linked to"* (<https://developer.chrome.com/blog/new-in-chrome-97>) — and `hidden=until-found`
  in 102, stating the problem it solves: *"content hidden in the collapsed sections becomes
  impossible to search using a find-in-page search"* and *"it isn't possible to link to text
  fragments inside collapsed regions"* (<https://developer.chrome.com/docs/css-ui/hidden-until-found>).

One more constraint on mixing routes: `:target` *"is set at document load and `history.back()`,
`history.forward()`, and `history.go()` method calls. But it is not changed when
`history.pushState()` and `history.replaceState()` methods are called"*
(<https://developer.mozilla.org/en-US/docs/Web/CSS/:target>), and `hashchange` *"does not fire when
the hash is modified using `history.pushState()` or `history.replaceState()`"*
(<https://developer.mozilla.org/en-US/docs/Web/API/Window/hashchange_event>). A `:target` scheme and
a `pushState` scheme cannot share one page.

### Mailing it

Gmail: *"For personal Gmail accounts, the limit is 25 MB"*, and over that *"Gmail automatically
removes the attachment and adds it as a Google Drive link"*
(<https://support.google.com/mail/answer/6584>). Its blocked-extension list runs from `.ade` through
`.img` and contains no `.htm` or `.html` entry
(<https://support.google.com/mail/answer/6590>). AP-1's page is 0.36 MB, so the cap is not near.

**No source found** for what a mail client does to an HTML attachment's contents, for whether a
multi-file set survives as an attachment, or for what a relative `<a href>` resolves to when only
one file of a set is attached. Those are assumptions, not facts, and #399 should treat them so.
**No source found** for any documented maximum size of a single HTML file.

### Route by route, from the above

| | Route 1, several files | Route 2, shell + routing | Route 3, one file, hidden sections |
| --- | --- | --- | --- |
| Opens from disk | yes, measured | **no** — the content fetch throws | yes |
| Back button | ordinary navigation, measured | fragment or query only | fragment only |
| Deep link | each page's own `file:` URL | fragment or query only | fragment only |
| Ctrl+F across the whole protocol | **no** — one page at a time | n/a | **only** with `<details>`, `hidden=until-found` or `content-visibility: auto` |
| Print | one command per page | n/a | one command, print CSS must open what is closed |
| Mail | several files, no source on what survives | n/a | one attachment |

## 3. What published bench-protocol formats do

Four were asked. **All four put one protocol in one document**, and all four chain to another
protocol by naming it, never by nesting it.

### protocols.io

A protocol is one object with one ordered step list. Its API's protocol object carries
`description, guidelines, before_start, warning, materials_text, materials, units, stats,
documents, funders, steps`, and a step carries `id, guid, previous_id, previous_guid, cases,
step, section, section_color` (<https://apidoc.protocols.io/>). **Sections are not containers** —
`section` is a label on each step, chained by `previous_guid`. The step component types include
`gotostep`, which jumps within one protocol. There is **no component type named "sub-protocol"**;
a `protocols` entity exists whose payload is a protocol object, and the docs do not say what it
renders as — **no source found**. Items carry a `fork_id` that *"identifies a parent protocol when
one exists"*, which is derivation, not nesting. Grouping is a separate type: the file-manager
filter takes *"protocol, collection, or document"*, type 1, 3 and 4.

A published export, measured: `picogram-input-multimodal-sequencing-pimms-cfaitice.pdf`
(<https://www.protocols.io/view/picogram-input-multimodal-sequencing-pimms-cfaitice.pdf>,
downloaded 2026-10-07) is **33 US-Letter pages in one continuous document**, produced by headless
Chrome — its PDF Producer is `Skia/PDF m149` and its Creator a `HeadlessChrome/149.0.0.0` string.
It has **no PDF outline, no named destinations and no table of contents**. Every page carries a
running footer: `protocols.io | https://dx.doi.org/10.17504/protocols.io.rm7vzywy5lx1/v1
January 25, 2023   N/33`. The DOI is version-pinned.

So the largest public protocol platform prints a protocol as one scroll, and gives paper no
navigation but "page N of 33" and an identifier.

### Nature Protocols

One article, with a fixed section sequence mandated in prose: *"The 'Introduction', 'Materials',
'Procedure', 'Timing', 'Troubleshooting' and 'Anticipated results' section headings are all
denoted with H1"* (<https://www.nature.com/nprot/for-authors/protocols>).

The numbering rule is the one worth carrying to #399: *"This must be provided as a numbered list
of direct experimental instructions"*, and *"If the protocol naturally breaks into separate
stages, then include subheadings and resume the numbered list."* **Step numbers run continuously
across sections; they do not reset.** Alternatives get a second alphabet — *"Letters of the Latin
alphabet (A), B), C)…) should be used to identify the different options, and Roman numerals i),
ii), iii)… should be used to break down the appropriate steps."*

Two devices key back into those numbers. Timing: *"include a timeline indicating the approximate
time each step or stage will take (e.g. Steps 1–3, 30 min; Steps 4+5, 2 h). Please also provide
this information as a list at the end of the procedure."* Troubleshooting: *"a table with the
column headings 'step', 'problem', 'possible reason', and 'solution'. The step number should be
given for where the problem is first observed."*

Chaining is citation, and formally a second article: a **Protocol Update** *"replaces the
associated Protocol"*; a **Protocol Extension** *"complements the associated Protocol (does not
replace it)"* and *"does not necessarily have authors in common"*
(<https://www.nature.com/nprot/for-authors/preparing-your-submission>). Authors are told to
reference others inline — *"You are encouraged, where appropriate, to reference other protocols,
including those published in Nature Protocols."*

The same page anticipates paper: *"Whilst Nature Protocols is an online product, some users will
print your protocol prior to use. Thus figures should be sized to be legible to users and to
facilitate printing."*

Whether the HTML article is one page or paginated was **not established** — no article was
fetched.

### Benchling

Its own words: *"An entry is the core document type in Benchling's Electronic Lab Notebook (ELN)"*
and *"Entries organize content into sections such as text, tables, file attachments, and Results
tables"* (<https://docs.benchling.com/docs/working-with-entries>). The API reads one entry as an
ordered part list — `DocumentPartSection`, `DocumentPartHeader`, `DocumentPartText` — so an entry
is one document of ordered parts. Reuse is a template, not a link: the warehouse carries
`entry_template`, *"Lists all notebook entry templates"*
(<https://docs.benchling.com/docs/wh-notebook>).

Its structured model is the interesting one for #395, because it is composition rather than
prose. In the Procedures beta: *"Procedure — A template for the recipe or assay"*, *"Method —
Unit operation"*, *"Parameter — Material input, material output, equipment, or parameter"*, and
*"Procedure run plan […] essentially a draft procedure run. It becomes a procedure run […] once
promoted"* (<https://docs.benchling.com/docs/wh-procedures>). There is an explicit flowchart of
edges *"determining which condition replicates flow into which condition replicates. No cycles."*
A procedure is built from unit operations with declared inputs and outputs, and template,
plan and run are three different things.

Benchling's user-facing help centre returned HTTP 403 to every fetch, so **no source found** for
what Benchling calls a protocol in the UI, whether protocols link, or how PDF export behaves.

### A vendor manual

Takara Bio's *In-Fusion Snap Assembly User Manual*, rev. 060822
(<https://www.takarabio.com/documents/User%20Manual/In/In-Fusion%20Snap%20Assembly%20User%20Manual.pdf>,
downloaded and measured 2026-10-07), is **one continuous 15-page PDF**, authored in Word.

Page 2 is a table of contents with page numbers and dotted leaders, followed by a Table of
Figures and a Table of Tables, each with page numbers. Sections are Roman numerals with lettered
sub-sections: I Introduction, II List of Components, III Additional Materials Required, IV PCR
Fragment Amplification, V Protocol, VI Transformation Procedure, VII Expected Results, VIII
Troubleshooting Guide, Appendix A, Appendix B. Every page footer reads `Page N of 15` plus the
revision stamp.

It navigates three ways, and each is worth naming:

1. **An outline table pointing at pages.** *"The table below is a general outline of the protocol
   […] Please refer to the specified pages for details on performing each step."*
2. **Cross-references by section number, inline.** *"Continue to the Transformation Procedure
   (Section VI)"*, *"see Section V.A for details"*, *"please refer to Appendix A"*.
3. **Prose pointers to a different document, with a catalogue number.** *"treat your PCR product
   with Cloning Enhancer (Cat. # 639615, Protocol not described in this document. Please refer to
   the Cloning Enhancer user manual)"*, and *"please read the In-Fusion Snap Assembly
   Multiple-Insert Cloning Protocol-At-A-Glance"*.

So the product has a **family** of documents — a full manual, a one-page "Protocol-At-A-Glance",
a variant manual, a web FAQ — not one document with a nav bar. Each is dated in its own filename
and footer.

NEB's E2621 manual and its web protocol pages, and Thermo Fisher's Gateway manual, all returned
HTTP 403 to every fetch attempted. **Not established** for those two vendors.

### What the four have in common

| | One document? | Chains by | Step numbers | Navigation on paper |
| --- | --- | --- | --- | --- |
| protocols.io | yes, flat steps with section labels | `gotostep` within; `fork_id`; Collections group | not established | none — no TOC, no bookmarks; `N/33` and a DOI per page |
| Nature Protocols | yes, fixed H1 sequence | citation; Update replaces, Extension complements | continuous, never reset; A)/i) for branches | Timing list and Troubleshooting table keyed by step number |
| Benchling | yes, one entry of ordered parts | composition: Procedure of Methods, no cycles; templates | not established | not established |
| Takara manual | yes, 15 pages | prose pointer naming the other manual and its catalogue number | per section (Roman + letter) | TOC with page numbers, outline table, section cross-refs, `Page N of 15` |

**Nobody found splits one protocol across pages.** What they split is one *project* into several
documents, each whole, each named, each citable — and the join is a sentence that names the other
document, not a link that needs a browser.

No format found uses breadcrumbs. No format found resets step numbers per section except the
vendor manual, which numbers sections instead of running one list.

## 4. The two-column layout #393 names

<https://zensical.org/docs/authoring/markdown/>, read and measured 2026-10-07.

### What it is, from the page itself

`<meta name="generator">` on that page reads **`zensical-0.0.68`**. The page loads exactly one
script, `../../assets/javascripts/bundle.3c842f4c.min.js`, and four stylesheets, of which three
are relative (`../../assets/stylesheets/modern/main.595f2bf5.min.css`, `palette…css`,
`extra.css`) and one is remote: `https://fonts.googleapis.com/css?family=Inter…`. Its inline
`#__config` block sets `"base": "../.."`. So the output is a directory of HTML files addressing
each other by relative path — route 1 — with one network dependency for fonts.

That same config lists the theme features in use:

```text
announce.dismiss, content.action.copy, content.action.edit, content.code.annotate,
content.code.copy, content.code.select, content.footnote.tooltips, content.tabs.link,
content.tooltips, navigation.footer, navigation.indexes, navigation.path,
navigation.sections, navigation.tabs, navigation.top, search.highlight
```

**`navigation.instant` is not in that list.** The reference layout does not use client-side
routing: every click on the left column is a full document load.

The config also points search at
`../../assets/javascripts/workers/search.7d14d953.min.js` — search runs in a Web Worker. Section
2 measured `new Worker` throwing under `file://`.

### How the two columns behave at width

Measured by emulating the viewport and reading each element's box:

| Viewport | Left column (`.md-sidebar--primary`) | Right column (`.md-sidebar--secondary`) | Content |
| --- | --- | --- | --- |
| 1426 px | 242 px wide, visible at x = 103 | 242 px wide, visible at x = 1081 | 736 px |
| 1100 px | 242 px wide at **x = −242**, `position: fixed` — off-screen, behind the hamburger | 242 px wide, visible at x = 858 | 826 px |
| 820 px, touch | off-screen, x = −242 | collapsed to a **38 × 38** box | 788 px |
| 390 px, mobile | off-screen, x = −242 | collapsed to 38 × 38 | 358 px, no horizontal page scroll |

The breakpoints come from the shipped stylesheet itself. `main.595f2bf5.min.css` is 148,595
bytes and holds 15 distinct media queries; the four that carry the layout are

| Query | Rule blocks | At a 16 px initial font |
| --- | --- | --- |
| `screen and (min-width:76.25em)` | 10 | 1220 px |
| `screen and (max-width:76.234375em)` | 9 | below 1220 px |
| `screen and (max-width:59.984375em)` | 6 | below 960 px |
| `screen and (max-width:44.984375em)` | 7 | below 720 px |

So the left column — the one that switches documents — is the first thing to go, and it goes at
1220 px. A 13-inch laptop at its default scaling, and every tablet, gets the drawer.

### What it does on paper

That stylesheet holds 20 `@media print` blocks. Four of them are the navigation:

```css
@media print{.md-header{display:none}}
@media print{.md-sidebar{display:none}}
@media print{.md-tabs{display:none}}
@media print{.md-feedback,.md-footer{display:none}}
```

Printing a page of that layout prints that page's content and nothing else — both columns
vanish, and with them every route to the other documents.

### What the projects say about themselves

Zensical is a static site generator by the Material for MkDocs authors: *"A modern static site
generator built by the creators of Material for MkDocs"*
(<https://github.com/zensical/zensical>, README). It *"shares the same core design principles
and philosophy"* (<https://zensical.org/docs/get-started/>) and it reads MkDocs configuration
(<https://zensical.org/docs/compatibility/mkdocs/>). It builds into a directory: *"This will
generate the static site in the configured site_dir, with the default being site"*
(<https://zensical.org/docs/usage/build/>). Both projects use the same class names, so Material's
documented mechanics apply, with the deltas noted below.

The breakpoints are declared in Material's source,
`src/templates/assets/stylesheets/_config.scss`:

```scss
$break-devices: (
  mobile:  (portrait: px2em(220px)  px2em(479.75px),  landscape: px2em(480px) px2em(719.75px)),
  tablet:  (portrait: px2em(720px)  px2em(959.75px),  landscape: px2em(960px) px2em(1219.75px)),
  screen:  (small:    px2em(1220px) px2em(1599.75px), medium: px2em(1600px) px2em(1999.75px),
            large:    px2em(2000px))
);
```

with `px2em($size, $base: 16px)`, which is where 1220 px, 960 px and 720 px come from. The theme
documents the first number in prose — top-level sections render in the sidebar *"for viewports
above `1220px`, but remain as-is on mobile"*
(<https://squidfunk.github.io/mkdocs-material/setup/setting-up-navigation/>). The 960 px boundary
is **in the source only**; no prose source found.

The columns have names: `_sidebar.scss` reads `// Primary sidebar with navigation` and
`// Secondary sidebar with table of contents`, and the same file carries the comments
`// [tablet -]: Show navigation as drawer` and `// [screen +]: Show navigation as sidebar`. The
drawer is opened by a CSS checkbox, `#__drawer`, with a `<label for="__drawer">` hamburger — no
script. Between 960 px and 1220 px, Material folds the table of contents into the drawer as a
second layer; Zensical instead turns it into a floating button at `right:.8rem;bottom:1.6rem`
opening a `max-height:50vh` popover. Above 960 px that button is hidden.

So the shape #393 names — left column switching documents, right column jumping within one —
**exists only above 1220 px**. From 960 to 1220 px it is one column and a hamburger. Below 960 px
it is neither.

### Whether the layout works from `file://`

It is the one thing both projects warn about, and in the same words. Material's offline plugin
page says:

> After building your project, switch to the `site` directory and open `index.html` in your
> browser […] However, you'll realize that the site search is gone. […] all calls to the Fetch
> API will error with a message like: `Cross origin requests are only supported for protocol
> schemes: http, [...]`

and that the plugin *"makes sure that site search keeps working by moving the search index to a
JavaScript file, and leveraging @squidfunk's `iframe-worker` shim"* and *"automatically disables
the `use_directory_urls` setting"* (<https://squidfunk.github.io/mkdocs-material/plugins/offline/>).
Zensical's page repeats it: *"all features that use the `fetch` API will error. Thus, when
building for offline usage, make sure to disable the following configuration settings: instant
navigation, site analytics, git repository, and comment systems"*
(<https://zensical.org/docs/setup/offline/>). Material's list adds versioning
(<https://squidfunk.github.io/mkdocs-material/setup/building-for-offline-usage/>).

MkDocs has a setting that exists for this case alone:

> This setting is needed when the documentation is hosted on systems that can't access the file
> `X/index.html` when given the URL `X`. […] For example, this needs to be set to `false` when:
> opening pages directly from the file system […]

(<https://www.mkdocs.org/user-guide/configuration/#use_directory_urls>; the default is `true`.)
Zensical sets it for you: *"this is automatically set to `false` when building for offline
usage"* (<https://zensical.org/docs/setup/basics/#use_directory_urls>).

Instant navigation is the client-side-routing route, and it is on both broken lists: *"clicks on
all internal links will be intercepted and dispatched via XHR without fully reloading the page
[…] Material for MkDocs now behaves like a Single Page Application"*, and *"you must set
`site_url` when using instant navigation, as instant navigation relies on the generated
`sitemap.xml`"* (<https://squidfunk.github.io/mkdocs-material/setup/setting-up-navigation/>).

One cost is easy to miss. Zensical's offline mode fetches its own shim from the network unless
you vendor it: *"This Javascript asset will be fetched from unpkg.com, unless you include it as
part of your own assets, in `extra.polyfills`"*, and *"The file name must contain the
`iframe-worker` substring, otherwise Zensical will fetch it again from unpkg.com"*
(<https://zensical.org/docs/setup/offline/>). The measured page also pulls fonts from
`fonts.googleapis.com`.

## 5. Printing a multi-page protocol

Printing a documentation site prints **one page per print**. Material's `_sidebar.scss` carries
`// [print]: Hide sidebar` with `@media print { display: none; }` on `.md-sidebar`, which covers
both columns; its compiled bundle holds 25 `@media print` blocks hiding the header, footer, tabs,
banner, dialog, clipboard buttons, feedback widget and header links, flattening content tabs to
`display: contents` and rendering annotations inline as `content: attr(data-md-annotation-id)`.
Zensical's bundle holds 20 such blocks, measured above. Neither project documents printing a
whole site; no source found.

Combining a site into one printable document is a third-party plugin's job.
`mkdocs-print-site-plugin` is *"MkDocs plugin that adds a print page to your site that combines
the entire site, allowing for easy export to PDF and standalone HTML"*, reached at
`/print_page/` or `print_page.html`, with a cover page and pagination options, and the
instruction *"Make sure to put `print-site` to the bottom of the plugin list"*
(<https://github.com/timvink/mkdocs-print-site-plugin>). It is not Material's and not Zensical's;
no source found that it runs under Zensical.

Read against section 1: splitting a protocol into pages means a bencher who wants paper either
prints each page separately and staples them, or the generator grows a second rendering path
whose only job is to put them back together. AP-1 is 66 printed pages today as one document.

### What the browser gives a print, in either shape

The page breaks in section 1's stylesheet (`break-inside: avoid` on figures, tables and the gel;
`break-after: avoid` on headings) are CSS Fragmentation properties and apply the same way in one
file or many. What differs between the shapes is only how many print commands the bencher issues
and whether anything renumbers across them. No source found for a browser printing several local
files in one command.

## 6. Open gaps

- **No source found** for Zensical's search index filename, or whether it is fetched or inlined;
  its docs say only that it is a new engine that *"doesn't use the Lunr pipeline"*
  (<https://zensical.org/docs/setup/search/>).
- **No source found** that `mkdocs-print-site-plugin` works under Zensical.
- **No source found** on whether Zensical has Material's `privacy` plugin, which *"will
  automatically download all external assets to distribute them with your documentation"*
  (<https://squidfunk.github.io/mkdocs-material/plugins/privacy/>) — relevant because the
  measured page loads fonts and the offline shim from the network.
- Material for MkDocs' maintenance status relative to Zensical appeared only in a secondary
  summary; **treat as unverified** until the announcement is read first-hand.
- The `file://` probe in section 2 is Chrome 154 only. Firefox and Safari were not measured, and
  Safari on iPadOS — a plausible bench device — was not measured at all. Everything found about
  opening local HTML on iPadOS was a forum post, so **no source found**: if a bench tablet
  matters to #399, it needs a real iPad, not a citation.
- **No source found** that Firefox and Safari auto-expand a closed `<details>` or support
  `hidden="until-found"` for find-in-page. Chrome's two announcements are sourced above; the
  other two engines were not confirmed first-hand.
- **No source found** for a verbatim statement that find-in-page skips plain `display: none`
  text. The two quotes above bound it — `content-visibility: hidden` *"must not be accessible"*
  and is *"similar to giving the contents `display: none`"*, and Chrome calls collapsed content
  *"impossible to search"* — but neither says it outright.
- **No source found** for a direct statement that hidden sections are omitted from print output.
  The nearest primary text is CSS Display on `display: none` — *"which causes the element's
  entire subtree to be left out of the box tree"* (<https://drafts.csswg.org/css-display/>) — no
  boxes, nothing to paint. Section 1 shows this repository already works around it by hand.
- **No source found** for protocols.io's "insert another protocol" behaviour, for whether its
  step numbers reset per section, or for whether a Nature Protocols HTML article is paginated.
  The protocols.io help centre is client-rendered and its API refuses anonymous reads.
- **Not established** for Benchling's user-facing protocol object or PDF export, and for NEB's
  and Thermo Fisher's manuals: all returned HTTP 403. A browser-driven fetch would settle them.
- The `data:` URI size ceiling has conflicting numbers in circulation; none was confirmed
  first-hand, and nothing here depends on it.
