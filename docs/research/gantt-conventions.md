---
search:
  exclude: true
---

# Gantt conventions: is there a ready way to chart a protocol whose durations are mostly unknown?

Research note for issue #404. Everything below was read or measured on **2026-10-07**.

Issue #404 measured one `liulab_mbio` bench protocol: of 42 steps, 17 hold a bounded duration
(70.0 h in total) and 24 hold none; four unattended overnight outgrowths are 74% of the held
hours; and no dependency, parallelism or hands-on time is recorded anywhere. The strings
"meanwhile", "in parallel", "next day" and "stopping point" appear zero times.

We decided not to draw a Gantt chart from that. This note asks whether we were wrong — whether
an established scheduling tool or convention already solves a chart in which many tasks have an
unknown duration and the dependency graph is absent.

## 1. The verdict, first

**No ready convention exists, and the evidence supports a plain schedule table with explicit
holes.**

Three findings decide it.

**Uncertainty collapses before anything is drawn.** Every established treatment of an unsure
duration — Microsoft Project's estimated-duration `?`, a PERT three-point estimate, a hammock
task — resolves to a single number before a bar reaches the page. The `?` is a flag on a value
that is still a value: Microsoft gives a new task "an estimate of one day" and draws the bar at
one day. PERT's `(O + 4M + P)/6` exists to collapse three numbers into one, and no mainstream
renderer draws the spread — only risk-analysis software does. Two products fade or lighten a bar
whose end is unknown, and neither is a convention anyone else follows (section 3). The
established answer to an unknown duration is to guess one, flag the guess, and draw the guess.

That is the one thing this repo has already refused to do. A `Hole` in
`liulab_mbio.protocol.model` is "a number nobody sourced, standing where the number would be",
and "filling one with a guess is the defect a hole exists to prevent". Adopting the
estimated-duration convention means adopting the defect.

**A Gantt with no predecessors is a bar chart with a date axis.** Nothing in a Gantt survives
the loss of the dependency graph except the bar lengths. With every task unconstrained and
starting as soon as possible, the bars stack at the project start; the chart's links, critical
path, float and sequence all read as empty. Section 3 finds no named convention for
"unsequenced" — the recognised answer is that the chart is not the right form yet.

**No renderer we could reach meets the output constraint and the data at once.** A chart here
must be inline SVG inside an offline `file://` page with no network and no script. Of the five
checked, only Vega-Lite through `vl-convert-python`, which we already ship, clears the output
constraint in process with no browser (section 4, measured). It can draw an uncertainty range
and a split bar. It cannot draw an open end: a null `x2` yields a rect with no width, so nothing
is drawn. An open end would be a convention we invent, not one we adopt.

**The venues that publish protocols chart none of this — but one of them has the field we are
missing.** Nature Protocols, protocols.io, Bio-protocol and JoVE all express a time plan as prose
or a list attached to a step or a step range, never as a schedule diagram (section 5). Nature's
form is `Steps 1-17, designing and cloning a targeted screen: 3-5 w, 1 w hands-on` — a duration
and a hands-on share, carried on a range of steps. That is a list, and it is the convention the
field settled on.

**What to build instead.** A schedule table: one row per step, its duration where the protocol
holds one, a `Hole` where it does not, and the kind of time it is. Filling in durations,
hands-on time and the first dependency is the work; the chart is a question to reopen once the
table stops being mostly holes.

## 2. How this was read and measured

- **Measured locally**: `vl_convert` 1.9.0 in this repo's pixi environment, driven from
  `pixi run python`. Four Gantt-shaped Vega-Lite specs were compiled with `vegalite_to_svg` and
  one inspected through `vegalite_to_scenegraph`. Every count in section 4 is from that run.
- **Read**: the `Step` and `Hole` dataclasses in `src/liulab_mbio/protocol/model.py`.
- **Read with a link**: tool documentation, library source and journal author instructions, as
  cited. Where a source is silent rather than negative, this note says so.
- The four non-Vega-Lite renderers were read, not run.

## 3. Uncertainty in a duration, and a plan with no dependencies

### What a tool does with a duration it does not know

| Device | What it is | Does it reach the bar? |
| --- | --- | --- |
| Estimated duration, the `?` suffix | A flag meaning "this number is a guess". Microsoft: the Estimated field "indicates whether the task's duration is flagged as an estimate", and "by default, when you first add a new task, it is given an estimate of one day" | No. The bar is drawn at the guessed value. The `?` shows in the Duration column |
| Elapsed duration, `ed` | Time counted around the clock: "7d is seven working days", "7ed is seven elapsed days, such as Monday - Sunday" | It changes the length, not the look |
| Manual scheduling and placeholder tasks | A task may hold a text duration "Project doesn't use for scheduling". A placeholder with only a duration gets "a lighter-colored placeholder Gantt bar"; with only one date it gets a "bookend" mark | **Yes** — the one mainstream tool that draws "not known" |
| PERT three-point estimate | Optimistic, most likely, pessimistic, collapsed by `(O + 4M + P)/6` | No. It is arithmetic before the drawing. Microsoft removed the PERT Analysis add-in in Project 2010 |
| Hammock task | A grouping whose duration is derived from the earliest start and latest finish of what it spans, never entered | No. It still resolves to one span |

Sources: Microsoft's
[Estimated field](https://support.microsoft.com/en-us/office/estimated-task-field-d5bd2410-9e1e-428f-bdb3-806cd2d8beff)
("a question mark appears with any estimated duration in the Duration field"), the
[Duration field](https://support.microsoft.com/en-us/project/duration-task-field), the
[DurationFormat reference](https://learn.microsoft.com/en-us/office-project/xml-data-interchange/durationformat-element)
for the elapsed units `em`, `eh`, `ed`, `ew`, `emo` and the estimated twins `d?` and `ed?`,
[how Project schedules tasks](https://support.microsoft.com/en-us/project/how-project-schedules-tasks-behind-the-scenes),
the [Placeholder field](https://support.microsoft.com/en-us/project/placeholder-task-field),
[changes in Project 2010](https://learn.microsoft.com/en-us/previous-versions/office/office-2010/cc178965(v=office.14))
listing "Pert Analysis" among removed add-ins, and
[hammock activity](https://en.wikipedia.org/wiki/Hammock_activity).

**Who draws a spread.** Only risk-analysis software. Intaver's RiskyProject shows "the ranges and
distributions for early and late start time and finish times" as "small triangles at the beginning
and end of each bar" ([Intaver](https://intaver.com/?p=3022)). That is a Monte Carlo result, which
needs a distribution per task and a dependency graph to propagate through. We have neither.

**Who draws an open end.** Two products, independently. Microsoft's lighter placeholder bar, above;
and BigPicture, where "if the start/end date is not defined, the relevant side of the taskbar is
faded"
([Appfire](https://appfire.atlassian.net/wiki/spaces/DLP/pages/2212497400)). OmniPlan answered a
2017 request for exactly this with "doesn't currently offer built-in functionality for
uncertain/estimated date entry"
([OmniGroup](https://discourse.omnigroup.com/t/uncertain-start-end-dates-for-a-task-way-to-show-it/34357)).
Two vendors with different pictures and no shared vocabulary is not a convention to adopt.

### A Gantt chart with no predecessors

**Every bar stacks at the start date.** As Soon As Possible is "the default in a project scheduled
from the start date"
([constraint type](https://support.microsoft.com/en-US/project/constraint-type-task-field)), and
"when you add a new task to a schedule, it automatically is scheduled to start on the project's
start date"; with nothing linked, "the project's duration is the same as the duration of the
longest task"
([scheduling](https://support.microsoft.com/en-us/project/how-project-schedules-tasks-behind-the-scenes)).
So 42 bars would all begin at hour zero, and the longest one would set the width.

**There is no name for it, because the field treats it as a defect, not a mode.** Schedule-quality
practice calls a task with no predecessor and no successor "dangling" or "missing logic", and the
DCMA 14-point assessment caps those at 5% of incomplete tasks
([Deltek](https://www.deltek.com/en/project-and-portfolio-management/project-scheduling/dcma-14-point-assessment),
[Tensix](https://tensix.com/the-14-point-assessment-and-schedule-missing-logic-inspection/)).
Our protocol would be at 100%. Float and the critical path, the two things a Gantt computes, are
undefined without logic.

**What is used instead is a list.** No standards body says so; the statement comes from vendors,
and the plainest is Instagantt's: a Gantt chart "without dependencies is just a list of tasks with
dates" ([Instagantt](https://www.instagantt.com/guides/project-dependencies)). The alternatives
offered are a task list, a calendar view or a Kanban board - each of which is a table with a time
column, drawn differently.

## 4. What the renderers can actually express

The constraint: inline SVG in a self-contained offline `file://` page, no network, no script at
view time.

| Renderer | (a) open end | (b) uncertainty range | (c) hands-on / unattended split | (d) static SVG, no JS |
| --- | --- | --- | --- | --- |
| **Vega-Lite** via `vl-convert-python` | No. A null `x2` draws nothing | Yes, `errorbar` or a layered `rule` | Yes, a stacked bar on one row | **Yes, in process, measured** |
| Mermaid `gantt` | No | No | No | Only through `mmdc`, which needs Node, Puppeteer and Chromium |
| frappe-gantt | No — a task with no end is dropped | No | Partial, one `progress` fill from the left edge | No render path outside a live browser DOM |
| vis-timeline | No — a missing end changes the shape to a box | No | No | No, and it draws HTML, not SVG |
| Plotly `px.timeline` | No, fakeable with a second faded bar | Yes, `error_x` | Yes, several bars on one y value | Only with Kaleido 1 and an installed Chrome |

**Vega-Lite, measured.** `vl_convert.vegalite_to_svg` compiled a layered Gantt-shaped spec in
0.50 s to 11,618 characters of SVG. The output starts with `<svg`, contains no `<script`, and
its only `http` strings are the SVG and xlink namespace URIs — nothing is fetched at view time.
A stacked hands-on and unattended split compiled, an `errorbar` layer compiled, and a dashed
`rule` with a `triangle-right` point compiled with `stroke-dasharray` intact. A `linearGradient`
survives, so a fading bar end is drawable.

**But Vega-Lite will not draw an unknown end for us.** With `x2` set to `null`, the scenegraph
holds the rect — `(x=75, x2=None, width=None)` — and it has no geometry, so nothing renders. The
caller must supply a finish. Any open end is a picture we design and a number we invent.

**One dependency note.** `vl-convert-python` is already a dependency, but
`src/liulab_mbio/plot/convert.py` uses only `svg_to_png` and `svg_to_pdf` — the package treats it
as an SVG rasteriser. `vegalite_to_svg` would be a new use of a shipped library, not a new
library.

Mermaid's `gantt` syntax has `section`, `milestone`, `excludes` and `until`, but a duration is a
single value such as `3d` and the only task tags are `active`, `done`, `crit` and `milestone`
([syntax](https://mermaid.js.org/syntax/gantt.html)); `mmdc` lists `puppeteer` as a peer
dependency ([mermaid-cli](https://github.com/mermaid-js/mermaid-cli)). frappe-gantt logs
`task "..." doesn't have an end date` and drops the task
([src/index.js](https://raw.githubusercontent.com/frappe/gantt/master/src/index.js)); its
`progress` field draws one `bar-progress` rect from the left
([src/bar.js](https://raw.githubusercontent.com/frappe/gantt/master/src/bar.js)). vis-timeline
"uses regular HTML DOM to render the timeline and items", and an item with no end becomes a box
rather than a range ([docs](https://visjs.github.io/vis-timeline/docs/timeline/)). Plotly marks
`x_start` and `x_end` required
([API](https://plotly.com/python-api-reference/generated/plotly.express.timeline.html)), carries
`error_x` on a bar trace ([reference](https://plotly.com/python/reference/bar/)), allows
"multiple bars on the same horizontal line" ([Gantt](https://plotly.com/python/gantt/)), and
exports SVG through Kaleido 1, which from v1 does not ship Chrome
([static export](https://plotly.com/python/static-image-export/)).

**Reading of the table.** Only one renderer clears the output constraint without a browser, and
even that one cannot say "unknown" — it can only draw a number we would have had to make up.

## 5. What the protocol venues publish

**None of the four draws a schedule.** All express time as prose or a list attached to a step or a
stage. None models a dependency, and only one names parallelism at all.

**Nature Protocols is the only venue with a real convention, and it is a list.** The author
guidelines ask authors to "include a timeline indicating the approximate time each step or stage
will take", to "also provide this information as a list at the end of the procedure", and, under
Procedure, to "include a TIMING callout with each subheading and state how long the section will
take to complete"
([for authors](https://www.nature.com/nprot/for-authors/protocols)). So the inline callout is per
stage, and step ranges appear only in the closing list.

The published form, read from Joung 2017 in
`reference_docs/synthesis_and_assembly/lentiviral-tolerance/assays/fm5-joung2017-natprotoc.txt`:

- inline, after a subheading: `Designing a custom sgRNA library ● TIMING 3-5 w; 1 w hands-on`
- in the closing Timing section: `Steps 1-17, designing and cloning a targeted screen: 3-5 w, 1 w hands-on`

Two things in that line matter to us. It carries **hands-on time as a second number on the same
duration** — the exact field issue #404 found missing — and it is attached to a **step range**,
not a step, which is how a protocol that cannot time every step still states a time.

The inline flags are five: `TIMING`, `PAUSE POINT` ("followed by a brief description of the
options available", for example "Can be left overnight at 4°C"), `CRITICAL STEP`, `CAUTION` and
`TROUBLESHOOTING`. A schedule diagram is not among the standard elements; the guidelines offer
only "consider including a flow diagram to demonstrate how the stages fit together".

**protocols.io stores a duration and nothing else about time.** A step's `components` array holds
an entity of type `duration`, whose data are `duration` (seconds), `label`, `pause`, `preset` and
`start` ([API](https://apidoc.protocols.io/)). Steps are an ordered list linked by `previous_guid`.
No dependency, no parallel step, no stopping-point field, and no protocol-level total.

**Bio-protocol asks for incubation times inside the sentence and nothing structured.** Its template
says to include "incubation time" per step and to "avoid using vague terms such as 'several,'
'enough,' 'about,' or 'around'", and it carries `Caution`, `Pause point`, `Critical` and
`See Troubleshooting` labels
([format](https://bio-protocol.org/bio101/standardformat.aspx)). There is no Timing section.

**JoVE records no timing at all.** Its author instructions give numbered steps of "1-2 actions" and
ask only: "Indicate any points at which the experiment can be paused and then restarted later"
([instructions](https://www.jove.com/files/Instructions_for_Authors.pdf)).

**The day-numbered schematic is a figure habit, not a standard.** Biology papers draw a line with
an anchor at day 0 and events at named days, and BioRender sells it as a template
([Mouse Experimental Timeline](https://www.biorender.com/template/mouse-experimental-timeline)). It
encodes where events sit on a calendar. It encodes no step duration, no hands-on split, no resource
and no dependency — which is to say it encodes none of what a Gantt chart is for.

## 6. What to build instead

A **schedule table**, one row per step, in the shape the venues already use:

| Column | Where it comes from |
| --- | --- |
| Step range | Nature's `Steps 1-17, ...: 3-5 w, 1 w hands-on` form — a range where no single step is timed |
| Duration | Summed from the step's `timers` and `programs` where it holds them |
| Hands-on | A new field. Nature publishes it, so it is a number a protocol is expected to know |
| Kind | Hands-on, attended wait, or unattended — the distinction that makes the four overnight outgrowths read as 74% of the hours without implying anyone is present for them |
| Missing | A `Hole`, where the protocol holds no duration |

A hole here is the same device the package already uses for an unsourced number, so 24 untimed
steps become 24 visible questions rather than 24 invented bars. `liulab_mbio.protocol.model.Step`
holds `timers`, `programs` and `holes` today and no duration, hands-on or predecessor field; the
table is what those fields would feed.

**Reopen the chart when the table is no longer mostly holes.** Two things would change the answer:
durations on most steps, and a first dependency — even "these two steps may run at once" as an
explicit pair. With both, Vega-Lite can draw it offline today, and the picture would be worth more
than the table. Without them, the chart would be an argument about invented numbers.
