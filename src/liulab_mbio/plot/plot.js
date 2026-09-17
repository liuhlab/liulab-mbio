// Hovering over anything a map draws shows its details, read from its group's data attributes;
// hovering over the notice of hidden labels lists them.
(() => {
  const tip = document.querySelector(".hover");
  const rows = ["name", "type", "span", "length"];

  const details = (item) =>
    item.dataset.hidden
      ? JSON.parse(item.dataset.hidden).map((label) => ["hidden", label])
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
      ? event.target.closest("[data-kind], [data-hidden]")
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
