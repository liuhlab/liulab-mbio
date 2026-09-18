// A drawing's page. Hovering over anything drawn shows its details, read from its group's data
// attributes, and hovering over a notice of hidden labels lists those that show. The switches
// show or hide each kind of item and each feature type in place, flip the map's shape, and show
// the sequence view. A click on an item highlights it in both views and scrolls the other to it.
// In the sequence view, hovering over a base shows its position, and a drag selects bases to copy,
// scrolling the view while the pointer lies past its top or bottom.
(() => {
  const tip = document.querySelector(".hover");
  const rows = ["name", "type", "span", "length"];
  const map = document.querySelector(".map");
  const view = document.querySelector(".sequence-view");

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
  };
  document.addEventListener("change", apply);
  // A browser may restore switches as they were left, rather than as the page was written.
  window.addEventListener("pageshow", apply);

  // The sequence view's rows: the bases each holds, counted as the stretch drawn counts them, and
  // how far down its strands and rail lie. Read when first needed.
  const whole = view?.querySelector(".rows");
  const length = Number(whole?.dataset.length);
  const cell = Number(whole?.dataset.cell);
  let found = null;
  const blocks = () =>
    (found ??= [...view.querySelectorAll(".row")].map((row) => ({
      row,
      start: Number(row.dataset.start),
      end: Number(row.dataset.end),
      top: Number(row.dataset.top),
      rail: Number(row.dataset.rail),
      bottom: Number(row.dataset.bottom),
    })));
  // Where a row's bases end: under the bottom strand, or at the rail when it hides.
  const foot = (block) => (view.classList.contains("one-strand") ? block.rail : block.bottom);
  // A position as a person reads it, from one counted along the stretch.
  const read = (at) => (at % length) + 1;

  // The base under a pointer, or with `near` the nearest base to it, or null off the strands.
  const baseAt = (event, near) => {
    const svg = view.querySelector("svg");
    const { x, y } = new DOMPoint(event.clientX, event.clientY).matrixTransform(
      svg.getScreenCTM().inverse(),
    );
    const all = blocks();
    let low = 0;
    let high = all.length - 1;
    while (low < high) {
      const middle = (low + high + 1) >> 1;
      if (all[middle].top <= y) low = middle;
      else high = middle - 1;
    }
    const block = all[low];
    const index = Math.floor(x / cell);
    const count = block.end - block.start;
    if (!near && (y < block.top || y > foot(block) || index < 0 || index >= count)) return null;
    return block.start + Math.min(Math.max(index, 0), count - 1);
  };

  // The selected bases, first and past the last, drawn under the rows and said in the caption.
  let chosen = null;
  let anchor = null;
  const marked = whole && document.createElementNS("http://www.w3.org/2000/svg", "path");
  if (marked) {
    marked.setAttribute("class", "chosen");
    whole.before(marked);
  }
  const overlapping = ([first, last]) =>
    blocks().filter((block) => block.start < last && first < block.end);
  const choose = (range) => {
    chosen = range;
    const parts = range
      ? overlapping(range).map((block) => {
          const left = (Math.max(range[0], block.start) - block.start) * cell;
          const right = (Math.min(range[1], block.end) - block.start) * cell;
          return `M${left} ${block.top}H${right}V${foot(block)}H${left}Z`;
        })
      : [];
    marked.setAttribute("d", parts.join(""));
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
    if (!view || !target?.closest(".plot svg")) return;
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
    view.querySelector("svg").setPointerCapture(event.pointerId);
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
  // A scroll under a still pointer, by the wheel or otherwise, moves the bases under it too.
  view?.addEventListener("scroll", () => {
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
    if (anchor !== null) event.preventDefault();
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
      scroll(map, there[0].getBoundingClientRect(), false);
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
    // A click in a caption, or one that selected bases, leaves the item highlighted.
    if (!target?.closest("svg") || onBase) return;
    const item = target.closest("[data-kind]");
    const chosenKey = item && key(item) !== selected ? key(item) : null;
    const groups = select(chosenKey);
    if (item && chosenKey) reveal(item, groups);
  });

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
