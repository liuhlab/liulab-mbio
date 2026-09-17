// A drawing's page. Hovering over anything drawn shows its details, read from its group's data
// attributes, and hovering over a notice of hidden labels lists those that show. The switches
// show or hide each kind of item and each feature type in place, and flip the map's shape. A click
// on an item highlights it wherever it is drawn.
(() => {
  const tip = document.querySelector(".hover");
  const rows = ["name", "type", "span", "length"];

  // The page's state lives in its switches: each input by its value.
  const switches = (name) =>
    new Map(
      [...document.querySelectorAll(`.switches input[name="${name}"]`)].map((input) => [
        input.value,
        input,
      ]),
    );
  const kinds = switches("kind");
  const types = switches("type");
  const shapes = switches("shape");

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
      for (const map of document.querySelectorAll(".map [data-shape]")) {
        map.hidden = map.dataset.shape !== shape;
      }
    }
  };
  document.querySelector(".switches")?.addEventListener("change", apply);
  // A browser may restore switches as they were left, rather than as the page was written.
  window.addEventListener("pageshow", apply);

  // An item is known by what hovering over it shows, wherever it is drawn.
  const key = (group) => JSON.stringify([group.dataset.kind, ...rows.map((row) => group.dataset[row])]);
  let selected = null;
  const select = (chosen) => {
    selected = chosen;
    for (const group of document.querySelectorAll(".plot [data-kind]")) {
      group.classList.toggle("selected", key(group) === chosen);
    }
  };
  document.querySelector(".plot").addEventListener("click", (event) => {
    const item = event.target instanceof Element ? event.target.closest("[data-kind]") : null;
    const chosen = item ? key(item) : null;
    select(chosen === selected ? null : chosen);
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

  let shown = null;
  document.addEventListener("pointermove", (event) => {
    const item = event.target instanceof Element
      ? event.target.closest(".plot [data-kind], .plot [data-hidden]")
      : null;
    if (!item) {
      tip.hidden = true;
      shown = null;
      return;
    }
    if (item !== shown) {
      tip.replaceChildren(
        ...details(item).map(([row, text]) => {
          const line = document.createElement("div");
          line.className = row;
          line.textContent = text;
          return line;
        }),
      );
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
