// A drawing's page. Hovering over anything drawn shows its details, read from its group's data
// attributes, and hovering over a notice of hidden labels lists those that show. The switches
// show or hide each kind of item and each feature type in place, flip the map's shape, and show
// the sequence view, beside the map or under it, in full rows or narrow ones, as the page's width
// allows. The map's caption zooms the shape shown, the circle freely and the line in steps, and a
// drag pans it. A click on an item highlights it in both views and scrolls the other to it. In the
// sequence view, hovering over a base shows its position, and a drag selects bases to copy,
// scrolling the view while the pointer lies past its top or bottom.
(() => {
  const tip = document.querySelector(".hover");
  const rows = ["name", "type", "span", "length"];
  const map = document.querySelector(".map");
  const view = document.querySelector(".sequence-view");
  // The sequence view's drawing at each row width, full rows first, and the one shown.
  const drawings = view ? [...view.querySelectorAll("[data-bases-per-row]")] : [];
  let drawn = drawings[0];
  const zoomer = map.querySelector("figcaption");
  const slider = zoomer?.querySelector('[name="zoom"]');

  // The map's shape shown, and whether it zooms: freely, or in steps between its drawings, each
  // laid out along a line twice as long as the last. Each shape keeps its zoom, in doublings, the
  // zoom a pinch last asked for, and where its pane was scrolled.
  const showing = () => map.querySelector(".shape:not([hidden])");
  const zooms = (shape) => Boolean(shape?.dataset.zoom);
  const steps = (shape) => [...shape.querySelectorAll(":scope > .step")];
  const image = (shape) =>
    (shape.querySelector(":scope > .step:not([hidden])") ?? shape).querySelector("svg");
  const kept = new Map();
  const held = (shape) => {
    if (!kept.has(shape)) kept.set(shape, { level: 0, aim: 0, left: 0, top: 0 });
    return kept.get(shape);
  };
  // Each step's drawing is shown at the scale its first is: as many times its first's width and
  // height as its view box is.
  const extent = (step) => step.querySelector("svg").viewBox.baseVal;
  for (const shape of map.querySelectorAll('[data-zoom="steps"]')) {
    const [first] = steps(shape).map(extent);
    for (const step of steps(shape)) {
      step.style.setProperty("--across", String(extent(step).width / first.width));
      step.style.setProperty("--down", String(extent(step).height / first.height));
    }
  }
  // A press on the map, which pans it once it moves more than a few pixels.
  let press = null;
  let panned = false;

  // The page's state lives in its switches: each input by its value.
  const switches = (name) =>
    new Map(
      [...document.querySelectorAll(`input[name="${name}"]`)].map((input) => [
        input.value,
        input,
      ]),
    );
  const kinds = switches("kind");
  const types = switches("type");
  const shapes = switches("shape");
  const sequence = switches("view").get("sequence");
  const strands = switches("strands").get("both");

  // Whether an item shows: its kind's switch is on, and a feature's type's too.
  const shows = ({ kind, type }) =>
    kinds.get(kind)?.checked !== false &&
    (kind !== "feature" || types.get(type)?.checked !== false);

  // The hidden labels a notice lists that show, and what it says of them, as `layers.notice` does.
  const listed = (notice) => JSON.parse(notice.dataset.hidden).filter(shows);
  const words = [
    ["cut_site", "enzyme site"],
    ["primer", "primer"],
    ["feature", "feature"],
  ];
  const said = (hidden) => {
    const parts = words
      .map(([kind, word]) => [hidden.filter((label) => label.kind === kind).length, word])
      .filter(([count]) => count)
      .map(([count, word]) => `${count} ${count === 1 ? word : `${word}s`}`);
    if (!parts.length) return "";
    const last = parts.pop();
    const all = parts.length ? `${parts.join(", ")} and ${last}` : last;
    return `${all} ${hidden.length === 1 ? "is" : "are"} hidden`;
  };

  const renotice = (notice) => {
    const text = notice.querySelector("text");
    const now = said(listed(notice));
    notice.classList.toggle("off", !now);
    if (!now || now === text.textContent) return;
    if (text.hasAttribute("textLength")) {
      // The words change length: keep their right edge where the map put it.
      const right = Number(text.getAttribute("x")) + Number(text.getAttribute("textLength"));
      text.setAttribute("x", String(right));
      text.setAttribute("text-anchor", "end");
      text.removeAttribute("textLength");
    }
    text.textContent = now;
  };

  const apply = () => {
    const was = showing();
    held(was).left = map.scrollLeft;
    held(was).top = map.scrollTop;
    for (const group of document.querySelectorAll(".plot [data-kind]")) {
      group.classList.toggle("off", !shows(group.dataset));
    }
    document.querySelectorAll(".plot .notice").forEach(renotice);
    const shape = [...shapes.values()].find((input) => input.checked)?.value;
    if (shape) {
      for (const one of map.querySelectorAll("[data-shape]")) {
        one.hidden = one.dataset.shape !== shape;
      }
    }
    if (view) {
      view.hidden = !sequence.checked;
      // The rows keep the room the bottom strand takes, so nothing moves.
      view.classList.toggle("one-strand", !strands.checked);
      if (chosen) choose(chosen);
    }
    arrange();
    const now = showing();
    if (zoomer) {
      zoomer.hidden = !zooms(now);
      slider.step = steps(now).length ? "1" : "any";
      slider.value = String(held(now).level);
    }
    if (now !== was) map.scrollTo(held(now).left, held(now).top);
    fit(now);
  };
  document.addEventListener("change", (event) => {
    if (event.target !== slider) apply();
  });
  // A browser may restore switches as they were left, rather than as the page was written.
  window.addEventListener("pageshow", apply);

  // Where the sequence view stands, and in which rows: the first of these arrangements that fits
  // the page's own width, given each drawing's natural width, its SVG's `width`.
  const plot = document.querySelector(".plot");
  const least = 320; // The narrowest the map goes beside the rows.
  const beside = (rows, width) => rows + width.gap + Math.min(width.map, least) <= width.room;
  const arrangements = [
    // Beside the map, full rows at their natural width, or else narrow ones.
    { under: false, narrow: false, fits: (width) => beside(width.full, width) },
    { under: false, narrow: true, fits: (width) => beside(width.narrow, width) },
    // Under the map, full rows at their natural width, or else narrow ones, scaled down where
    // they do not fit.
    { under: true, narrow: false, fits: (width) => width.full <= width.room },
    { under: true, narrow: true, fits: () => true },
  ];
  const natural = (drawing) => Number(drawing.querySelector("svg").getAttribute("width"));
  const arrange = () => {
    if (!view) return;
    const style = getComputedStyle(plot);
    const width = {
      room: plot.clientWidth - parseFloat(style.paddingLeft) - parseFloat(style.paddingRight),
      gap: parseFloat(style.columnGap),
      map: natural(map.querySelector("[data-shape]:not([hidden])")),
      full: natural(drawings[0]),
      narrow: natural(drawings.at(-1)),
    };
    const { under, narrow } = arrangements.find((one) => one.fits(width));
    plot.classList.toggle("under", under && !view.hidden);
    drawn = narrow ? drawings.at(-1) : drawings[0];
    for (const one of drawings) one.hidden = one !== drawn;
    if (!view.hidden) keep();
  };
  new ResizeObserver(arrange).observe(plot);

  // Each drawing's rows: the bases each holds, counted as the stretch drawn counts them, and how
  // far down its strands and rail lie. Read when first needed.
  const whole = view?.querySelector(".rows");
  const length = Number(whole?.dataset.length);
  const cell = Number(whole?.dataset.cell);
  const found = new Map();
  const blocks = (drawing = drawn) => {
    if (!found.has(drawing)) {
      found.set(
        drawing,
        [...drawing.querySelectorAll(".row")].map((row) => ({
          row,
          start: Number(row.dataset.start),
          end: Number(row.dataset.end),
          top: Number(row.dataset.top),
          rail: Number(row.dataset.rail),
          bottom: Number(row.dataset.bottom),
        })),
      );
    }
    return found.get(drawing);
  };
  // Where a row's bases end: under the bottom strand, or at the rail when it hides.
  const foot = (block) => (view.classList.contains("one-strand") ? block.rail : block.bottom);
  // A position as a person reads it, from one counted along the stretch.
  const read = (at) => (at % length) + 1;
  // How many of the rows shown, from the first, `before` holds for: it holds for every row before
  // one it does not hold for.
  const leading = (before) => {
    const all = blocks();
    let low = 0;
    let high = all.length;
    while (low < high) {
      const middle = (low + high) >> 1;
      if (before(all[middle])) low = middle + 1;
      else high = middle;
    }
    return low;
  };
  // The row shown at an index, or the nearest there is.
  const blockAt = (index) => blocks()[Math.min(Math.max(index, 0), blocks().length - 1)];
  // A client point in the rows shown, in the points they are laid out in.
  const inRows = (clientX, clientY) =>
    new DOMPoint(clientX, clientY).matrixTransform(
      drawn.querySelector("svg").getScreenCTM().inverse(),
    );

  // The base under a pointer, or with `near` the nearest base to it, or null off the strands.
  const baseAt = (event, near) => {
    const { x, y } = inRows(event.clientX, event.clientY);
    // The last row starting above the pointer, or the first.
    const block = blockAt(leading((one) => one.top <= y) - 1);
    const index = Math.floor(x / cell);
    const count = block.end - block.start;
    if (!near && (y < block.top || y > foot(block) || index < 0 || index >= count)) return null;
    return block.start + Math.min(Math.max(index, 0), count - 1);
  };

  // The selected bases, first and past the last, drawn under each drawing's rows and said in the
  // caption.
  let chosen = null;
  let anchor = null;
  const marks = new Map(
    drawings.map((drawing) => {
      const marked = document.createElementNS("http://www.w3.org/2000/svg", "path");
      marked.setAttribute("class", "chosen");
      drawing.querySelector(".rows").before(marked);
      return [drawing, marked];
    }),
  );
  const overlapping = ([first, last], drawing = drawn) =>
    blocks(drawing).filter((block) => block.start < last && first < block.end);
  const choose = (range) => {
    chosen = range;
    for (const [drawing, marked] of marks) {
      const parts = range
        ? overlapping(range, drawing).map((block) => {
            const left = (Math.max(range[0], block.start) - block.start) * cell;
            const right = (Math.min(range[1], block.end) - block.start) * cell;
            return `M${left} ${block.top}H${right}V${foot(block)}H${left}Z`;
          })
        : [];
      marked.setAttribute("d", parts.join(""));
    }
    const [first, last] = range ?? [0, 0];
    view.querySelector(".selection").textContent = range
      ? `${read(first)} .. ${read(last - 1)} (${last - first} bp)`
      : "";
    view.querySelector(".copy").hidden = !range;
  };
  // The top strand's bases selected.
  const bases = () =>
    overlapping(chosen)
      .map((block) =>
        block.row
          .querySelector(".strand.top text")
          .textContent.slice(
            Math.max(chosen[0], block.start) - block.start,
            Math.min(chosen[1], block.end) - block.start,
          ),
      )
      .join("");

  // The part of a view in sight: under its caption, which covers what scrolls beneath it.
  const sight = (figure) => {
    const frame = figure.getBoundingClientRect();
    return {
      top: frame.top + (figure.querySelector("figcaption")?.offsetHeight ?? 0),
      bottom: frame.top + figure.clientHeight,
      left: frame.left,
      right: frame.left + figure.clientWidth,
    };
  };

  // The first base of the top row in sight, as scrolling last left it, for a change of rows to keep
  // at the top.
  let heading = null;
  // The top row in sight: the first whose bases reach below the caption.
  const topmost = () => {
    const { left, top } = sight(view);
    const { y } = inRows(left, top);
    return blockAt(leading((block) => foot(block) <= y));
  };
  const holds = (block, at) => block.start <= at && at < block.end;
  // Scrolls the row holding that base to the top, unless the top row in sight holds it.
  const keep = () => {
    if (heading === null || holds(topmost(), heading)) return;
    const block = blockAt(leading((one) => one.end <= heading));
    if (block === blocks()[0]) view.scrollTop = 0;
    else scroll(view, block.row.getBoundingClientRect(), true);
  };

  // While a drag lasts: where the pointer last was, and the frame loop scrolling the view.
  let pointer = null;
  let loop = 0;
  // Extends the selection from its anchor to the nearest base in sight at the pointer.
  const extend = () => {
    const { top, bottom } = sight(view);
    const clientY = Math.min(Math.max(pointer.clientY, top), bottom);
    const at = baseAt({ clientX: pointer.clientX, clientY }, true);
    const range = [Math.min(anchor, at), Math.max(anchor, at) + 1];
    if (range[0] !== chosen?.[0] || range[1] !== chosen?.[1]) choose(range);
  };
  // Each second, a pointer past the rows in sight scrolls the view this many times as far.
  const pace = 20;
  const glide = (then, owed) => {
    loop = requestAnimationFrame((now) => {
      const { top, bottom } = sight(view);
      const y = pointer.clientY;
      const past = y < top ? y - top : Math.max(y - bottom, 0);
      // Whole pixels only, so a slow scroll is not lost to rounding.
      const due = past ? owed + (past * pace * (now - (then ?? now))) / 1000 : 0;
      const step = Math.trunc(due);
      if (step) {
        view.scrollTop += step;
        extend();
      }
      glide(now, due - step);
    });
  };

  let onBase = false;
  document.addEventListener("pointerdown", (event) => {
    const target = event.target instanceof Element ? event.target : null;
    onBase = false;
    if (!view || !target?.closest(".sequence-view svg")) return;
    const at =
      event.button === 0 && target.closest(".sequence-view") && !target.closest("[data-kind]")
        ? baseAt(event, false)
        : null;
    if (at === null) {
      if (chosen) choose(null);
      return;
    }
    onBase = true;
    anchor = at;
    pointer = event;
    drawn.querySelector("svg").setPointerCapture(event.pointerId);
    document.getSelection()?.removeAllRanges();
    choose([at, at + 1]);
    cancelAnimationFrame(loop);
    glide(null, 0);
  });
  document.addEventListener("pointermove", (event) => {
    if (anchor === null) return;
    pointer = event;
    extend();
  });
  view?.addEventListener("scroll", () => {
    const block = topmost();
    if (heading === null || !holds(block, heading)) heading = block.start;
    // And a scroll under a still pointer, by the wheel or otherwise, moves the bases under it too.
    if (anchor !== null) extend();
  });
  const release = () => {
    anchor = null;
    cancelAnimationFrame(loop);
  };
  document.addEventListener("pointerup", release);
  document.addEventListener("pointercancel", release);
  // A drag over the bases selects them, not the text they are set in.
  document.addEventListener("selectstart", (event) => {
    if (anchor !== null || press) event.preventDefault();
  });

  // Copying puts the selected bases on the clipboard, unless text is selected, but for the button.
  let copying = false;
  document.addEventListener("copy", (event) => {
    const text = document.getSelection()?.isCollapsed === false;
    if (!chosen || (text && !copying)) return;
    event.clipboardData.setData("text/plain", bases());
    event.preventDefault();
  });
  view?.querySelector(".copy").addEventListener("click", () => {
    copying = true;
    document.execCommand("copy");
    copying = false;
  });

  // An item is known by what hovering over it shows, wherever it is drawn.
  const key = (group) => JSON.stringify([group.dataset.kind, ...rows.map((row) => group.dataset[row])]);
  let selected = null;
  // Marks each group of the chosen item, and returns them.
  const select = (chosenKey) => {
    selected = chosenKey;
    const groups = [];
    for (const group of document.querySelectorAll(".plot [data-kind]")) {
      const mine = key(group) === chosenKey;
      group.classList.toggle("selected", mine);
      if (mine) groups.push(group);
    }
    return groups;
  };

  // Scrolls a view so `box` shows: at its top, or when `start` is false, centred if it does not.
  const scroll = (figure, box, start) => {
    const { top, bottom, left, right } = sight(figure);
    if (start || box.top < top || box.bottom > bottom) {
      const offset = start ? box.top - top : (box.top + box.bottom - top - bottom) / 2;
      figure.scrollTop += offset;
    }
    if (box.left < left || box.right > right) {
      figure.scrollLeft += (box.left + box.right - left - right) / 2;
    }
  };

  // Scrolls the view other than the one clicked to the item: the map to where it is drawn, and the
  // sequence view to the row holding its start, or the first row it lies in.
  const reveal = (clicked, groups) => {
    const other = clicked.closest(".map") ? view : map;
    if (!other || other.hidden) return;
    const there = groups.filter((group) => other.contains(group) && !group.closest("[hidden]"));
    if (!there.length) return;
    if (other === map) {
      // Its name, which a zoomed map may show when the whole of a long arrow cannot.
      const name = there.map((group) => group.querySelector("text")).find(Boolean);
      scroll(map, (name ?? there[0]).getBoundingClientRect(), false);
      return;
    }
    const at = parseInt(clicked.dataset.span, 10) - 1;
    const held = blocks().filter((block) => there.some((group) => block.row.contains(group)));
    const block =
      held.find(({ start, end }) => [at, at + length].some((one) => start <= one && one < end)) ??
      held[0];
    scroll(view, block.row.getBoundingClientRect(), true);
  };

  document.querySelector(".plot").addEventListener("click", (event) => {
    const target = event.target instanceof Element ? event.target : null;
    // A click in a caption, one that selected bases, or a pan, leaves the item highlighted.
    if (!target?.closest("svg") || onBase || panned) return;
    if (chosen && target.closest(".map")) choose(null);
    const item = target.closest("[data-kind]");
    const chosenKey = item && key(item) !== selected ? key(item) : null;
    const groups = select(chosenKey);
    if (item && chosenKey) reveal(item, groups);
  });

  // Where the shape shown was scrolled, kept as its pane scrolls.
  map.addEventListener("scroll", () => {
    held(showing()).left = map.scrollLeft;
    held(showing()).top = map.scrollTop;
  });

  // Measures a shape at 1×: the room its drawing takes, which it keeps while zoomed, and the room
  // its pane holds right of and below the drawing, which stays there. So the pane scrolls as much
  // further as the drawing grew, and any point of it can stay where it is.
  const pin = (shape) => {
    const box = shape.querySelector("svg").getBoundingClientRect();
    const frame = map.getBoundingClientRect();
    const right = box.right - frame.left - map.clientLeft + map.scrollLeft;
    const below = box.bottom - frame.top - map.clientTop + map.scrollTop;
    shape.style.setProperty("--width", `${box.width}px`);
    shape.style.setProperty("--height", `${box.height}px`);
    shape.style.setProperty("--right", `${Math.max(map.scrollWidth - right, 0)}px`);
    shape.style.setProperty("--below", `${Math.max(map.scrollHeight - below, 0)}px`);
    return box;
  };
  // Measures a zoomed shape again, at its first step when it steps, and scales where the pane was
  // scrolled by how much it changed.
  const fit = (shape) => {
    if (!shape.classList.contains("zoomed")) return;
    const by = 1 / parseFloat(shape.style.getPropertyValue("--width"));
    const [first, shown] = [steps(shape)[0], steps(shape)[held(shape).level]];
    shape.classList.remove("zoomed");
    if (shown) [shown.hidden, first.hidden] = [true, false];
    const { width } = pin(shape);
    if (shown) [first.hidden, shown.hidden] = [true, false];
    shape.classList.add("zoomed");
    map.scrollTo(held(shape).left * width * by, held(shape).top * width * by);
  };
  new ResizeObserver(() => fit(showing())).observe(plot);

  // Where the last zoom scrolled the pane, before the pane rounded it to whole pixels, so that
  // zooming step by step does not add up the rounding.
  let wanted = null;
  // Zooms the shape shown to `level` doublings, a whole step for a shape that steps, keeping the
  // point at (x, y) where it is. The drawing grows from its top left, or is swapped for the step's,
  // and the pane scrolls over it.
  const zoom = (level, x, y) => {
    const shape = showing();
    if (!zooms(shape)) return;
    const all = steps(shape);
    const most = all.length ? all.length - 1 : Number(slider.max);
    const to = Math.min(Math.max(all.length ? Math.round(level) : level, 0), most);
    const by = 2 ** (to - held(shape).level);
    const was = image(shape);
    const box = held(shape).level ? was.getBoundingClientRect() : pin(shape);
    const { scrollLeft, scrollTop } = map;
    const near = (one, other) => Math.abs(one - other) < 1;
    const [left, top] =
      wanted?.shape === shape && near(wanted.left, scrollLeft) && near(wanted.top, scrollTop)
        ? [wanted.left, wanted.top]
        : [scrollLeft, scrollTop];
    // How far the point lies from the drawing's top left, had the pane been scrolled to exactly
    // `left` and `top`, and how far it will.
    const at = [x - box.left - scrollLeft + left, y - box.top - scrollTop + top];
    let then = at.map((one) => one * by);
    if (all.length) {
      // Each step's line runs from x = 0 twice as far as the last's, its text the same size: the
      // point moves along with the line, and stays as high above it.
      const scale = parseFloat(shape.style.getPropertyValue("--width")) / extent(all[0]).width;
      const [from, into] = [was.viewBox.baseVal, extent(all[to])];
      const [u, v] = [at[0] / scale + from.x, at[1] / scale + from.y];
      then = [(u * by - into.x) * scale, (v - into.y) * scale];
      all.forEach((step, index) => {
        step.hidden = index !== to;
      });
    }
    shape.style.setProperty("--zoom", String(2 ** to));
    shape.classList.toggle("zoomed", to > 0);
    const [goal, fall] = [left + then[0] - at[0], top + then[1] - at[1]];
    map.scrollTo(goal, fall);
    // As far as the pane scrolls, so that an end it stopped at is not taken for rounding.
    const within = (one, room) => Math.min(Math.max(one, 0), room);
    wanted = {
      shape,
      left: within(goal, map.scrollWidth - map.clientWidth),
      top: within(fall, map.scrollHeight - map.clientHeight),
    };
    held(shape).level = held(shape).aim = to;
    slider.value = String(to);
  };

  // The slider zooms about the middle of the drawing in sight, and Ctrl or ⌘ with the wheel, as a
  // trackpad's pinch sends it, about the pointer. A plain wheel scrolls.
  slider?.addEventListener("input", () => {
    const { top, bottom, left, right } = sight(map);
    const box = image(showing()).getBoundingClientRect();
    const x = (Math.max(left, box.left) + Math.min(right, box.right)) / 2;
    const y = (Math.max(top, box.top) + Math.min(bottom, box.bottom)) / 2;
    zoom(Number(slider.value), x, y);
  });
  zoomer?.querySelector(".reset").addEventListener("click", () => {
    zoom(0, 0, 0);
    map.scrollTo(0, 0);
  });
  map.addEventListener(
    "wheel",
    (event) => {
      const shape = showing();
      if (!(event.ctrlKey || event.metaKey) || !zooms(shape)) return;
      event.preventDefault();
      // A pinch sends 100 × ln(its scale); a wheel's notch zooms at most half a doubling.
      const lines = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? map.clientHeight : 1;
      const by = Math.min(Math.max(-(event.deltaY * lines) / 100 / Math.LN2, -0.5), 0.5);
      const hold = held(shape);
      if (!steps(shape).length) {
        zoom(hold.level + by, event.clientX, event.clientY);
        return;
      }
      // On a shape that steps, a notch steps once. A pinch's many small steps add up, and the step
      // nearest their sum shows, a tie going the way they moved.
      const asked = Math.abs(by) === 0.5 ? hold.level + Math.sign(by) : hold.aim + by;
      const aim = Math.min(Math.max(asked, 0), Number(slider.max));
      zoom(by > 0 ? Math.floor(aim + 0.5) : Math.ceil(aim - 0.5), event.clientX, event.clientY);
      hold.aim = aim;
    },
    { passive: false },
  );

  // A press on the map that moves more than a few pixels pans it, and is not a click. A touch
  // pans as the browser scrolls.
  document.addEventListener("pointerdown", (event) => {
    const target = event.target instanceof Element ? event.target : null;
    const { clientX: x, clientY: y, pointerId: id } = event;
    panned = false;
    press =
      event.button === 0 && event.pointerType !== "touch" && target?.closest(".map svg")
        ? { x, y, id, left: map.scrollLeft, top: map.scrollTop }
        : null;
  });
  document.addEventListener("pointermove", (event) => {
    if (!press) return;
    const x = event.clientX - press.x;
    const y = event.clientY - press.y;
    if (!panned && Math.hypot(x, y) <= 4) return;
    if (!panned) {
      panned = true;
      map.setPointerCapture(press.id);
      map.classList.add("panning");
      tip.hidden = true;
    }
    map.scrollTo(press.left - x, press.top - y);
  });
  const drop = () => {
    press = null;
    map.classList.remove("panning");
  };
  document.addEventListener("pointerup", drop);
  document.addEventListener("pointercancel", drop);

  const details = (item) =>
    item.dataset.hidden
      ? listed(item).map((label) => ["hidden", label.label])
      : rows.filter((row) => item.dataset[row]).map((row) => [row, item.dataset[row]]);

  const place = (event) => {
    const gap = 14;
    const x = event.clientX + gap + tip.offsetWidth > window.innerWidth
      ? event.clientX - gap - tip.offsetWidth
      : event.clientX + gap;
    const y = event.clientY + gap + tip.offsetHeight > window.innerHeight
      ? event.clientY - gap - tip.offsetHeight
      : event.clientY + gap;
    tip.style.left = `${Math.max(0, x)}px`;
    tip.style.top = `${Math.max(0, y)}px`;
  };

  const say = (said) =>
    tip.replaceChildren(
      ...said.map(([row, text]) => {
        const line = document.createElement("div");
        line.className = row;
        line.textContent = text;
        return line;
      }),
    );

  let shown = null;
  document.addEventListener("pointermove", (event) => {
    if (press && panned) return;
    const target = event.target instanceof Element ? event.target : null;
    const item = target?.closest(".plot [data-kind], .plot [data-hidden]");
    const base = !item && target?.closest(".sequence-view svg") ? baseAt(event, false) : null;
    if (!item && base === null) {
      tip.hidden = true;
      shown = null;
      return;
    }
    if (base !== null) {
      say([["position", `position ${read(base)}`]]);
      shown = null;
    } else if (item !== shown) {
      say(details(item));
      shown = item;
    }
    tip.hidden = false;
    place(event);
  });
  document.addEventListener("pointerleave", () => {
    tip.hidden = true;
    shown = null;
  });
})();
