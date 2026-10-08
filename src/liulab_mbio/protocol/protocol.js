(function () {
  "use strict";

  var storageKey = "liulab-protocol:" + (document.body.getAttribute("data-protocol") || "");

  // Storage can be missing or throw (private windows, blocked site data); the page works without it.
  function load() {
    try {
      return JSON.parse(window.localStorage.getItem(storageKey) || "{}") || {};
    } catch (error) {
      return {};
    }
  }

  function save() {
    try {
      window.localStorage.setItem(storageKey, JSON.stringify(state));
    } catch (error) {
      /* checks last for this visit only */
    }
  }

  var state = load();

  // One step is a step: the count reads as a sentence, as render.py writes it.
  function steps(n) {
    return n + (n === 1 ? " step" : " steps");
  }

  function all(selector, root) {
    return Array.prototype.slice.call((root || document).querySelectorAll(selector));
  }

  // Check marks.
  var boxes = all('input[type="checkbox"][data-key]');
  var stepBoxes = all("input.done");
  var progress = document.querySelector(".progress");

  function refresh() {
    var done = 0;
    stepBoxes.forEach(function (box) {
      var step = box.closest(".step");
      if (step) step.classList.toggle("is-done", box.checked);
      if (box.checked) done += 1;
    });
    if (progress) progress.textContent = done + " of " + steps(stepBoxes.length) + " done";
  }

  boxes.forEach(function (box) {
    box.checked = state[box.getAttribute("data-key")] === true;
    box.addEventListener("change", function () {
      var key = box.getAttribute("data-key");
      if (box.checked) state[key] = true;
      else delete state[key];
      save();
      refresh();
    });
  });
  refresh();

  // A run's index: how far the bench got in each protocol, read from that page's own store.
  // Every file:// page shares one store, and each page keys by its own content, so the index
  // reads the marks without the protocol pages writing anything twice.
  all("[data-page-key]").forEach(function (item) {
    var total = Number(item.getAttribute("data-steps")) || 0;
    var label = item.querySelector(".page-progress");
    var marks;
    try {
      marks = JSON.parse(
        window.localStorage.getItem("liulab-protocol:" + item.getAttribute("data-page-key")) || "{}"
      );
    } catch (error) {
      return;
    }
    if (!marks || !label || !total) return;
    var done = 0;
    for (var n = 1; n <= total; n += 1) if (marks["step-" + n] === true) done += 1;
    label.textContent = done + " of " + steps(total) + " done";
    item.classList.toggle("is-started", done > 0);
  });

  // Everything this page remembers is one object, so resetting it is emptying that object and
  // reading the page again: marks, reaction counts and timers all go back to what was written.
  var reset = document.querySelector("button.reset");
  if (reset) {
    reset.addEventListener("click", function () {
      Object.keys(state).forEach(function (key) { delete state[key]; });
      save();
      window.location.reload();
    });
  }

  var print = document.querySelector("button.print");
  if (print) print.addEventListener("click", function () { window.print(); });

  // Reaction tables: the arithmetic render.py writes the mix column with, and the rule of
  // protocol.model.number, so the column reads the same after the count changes.
  var SUPERSCRIPT = "⁰¹²³⁴⁵⁶⁷⁸⁹";

  function number(value) {
    if (value === 0) return "0";
    var parts = value.toExponential(2).split("e");
    var power = Number(parts[1]);
    if (Math.abs(value) >= 0.001 && Math.abs(value) < 1e6) {
      var fixed = value.toFixed(Math.max(0, 2 - power));
      if (fixed.indexOf(".") >= 0) fixed = fixed.replace(/0+$/, "").replace(/\.$/, "");
      var halves = fixed.split(".");
      halves[0] = halves[0].replace(/\B(?=(\d{3})+(?!\d))/g, ",");
      return halves.join(".");
    }
    var mantissa = parts[0].replace(/\.?0+$/, "");
    var raised = String(power).replace(/\d/g, function (d) { return SUPERSCRIPT.charAt(d); });
    return mantissa + " × 10" + raised.replace("-", "⁻");
  }

  all(".reaction").forEach(function (figure) {
    var input = figure.querySelector("input.rxn-count");
    if (!input) return;
    var overage = parseFloat(figure.getAttribute("data-overage")) || 0;
    var key = input.getAttribute("data-key");

    function update() {
      var reactions = Math.max(1, Math.floor(Number(input.value) || 1));
      var scale = reactions * (1 + overage);
      all("[data-ul]", figure).forEach(function (cell) {
        cell.textContent = number(parseFloat(cell.getAttribute("data-ul")) * scale);
      });
      all(".rxn-n", figure).forEach(function (span) { span.textContent = String(reactions); });
      return reactions;
    }

    // A browser restores a typed-in value across a reload, so the count is set either way: what
    // the page remembers, or the count render.py wrote.
    if (typeof state[key] === "number") input.value = String(state[key]);
    else input.value = input.getAttribute("value") || "1";
    update();
    input.addEventListener("input", function () {
      state[key] = update();
      save();
    });
  });

  // Copy buttons.
  function fallbackCopy(text) {
    var area = document.createElement("textarea");
    area.value = text;
    area.setAttribute("readonly", "");
    area.style.position = "fixed";
    area.style.opacity = "0";
    document.body.appendChild(area);
    area.select();
    try {
      return document.execCommand("copy");
    } catch (error) {
      return false;
    } finally {
      document.body.removeChild(area);
    }
  }

  function flash(button, label) {
    var original = button.getAttribute("data-label") || button.textContent;
    button.setAttribute("data-label", original);
    button.textContent = label;
    button.classList.add("is-copied");
    window.setTimeout(function () {
      button.textContent = original;
      button.classList.remove("is-copied");
    }, 1400);
  }

  all("button.copy").forEach(function (button) {
    button.addEventListener("click", function () {
      var text = button.getAttribute("data-copy") || "";
      var done = function () { flash(button, "Copied"); };
      var failed = function () { flash(button, fallbackCopy(text) ? "Copied" : "Copy failed"); };
      if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(text).then(done, failed);
      } else {
        failed();
      }
    });
  });

  // Timers.
  function clock(seconds) {
    var s = Math.max(0, Math.round(seconds));
    var h = Math.floor(s / 3600);
    var m = Math.floor((s % 3600) / 60);
    var rest = String(s % 60).padStart(2, "0");
    return h ? h + ":" + String(m).padStart(2, "0") + ":" + rest : m + ":" + rest;
  }

  function alarm() {
    try {
      var Context = window.AudioContext || window.webkitAudioContext;
      var context = new Context();
      [0, 0.35, 0.7].forEach(function (offset) {
        var tone = context.createOscillator();
        var gain = context.createGain();
        tone.frequency.value = 880;
        gain.gain.value = 0.15;
        tone.connect(gain);
        gain.connect(context.destination);
        tone.start(context.currentTime + offset);
        tone.stop(context.currentTime + offset + 0.2);
      });
    } catch (error) {
      /* no sound available */
    }
    if (navigator.vibrate) navigator.vibrate([200, 100, 200]);
  }

  // A running timer keeps its deadline and a paused one the seconds it has left, so turning the
  // page or closing the tab does not lose an incubation. A deadline already past comes back
  // finished and silent: the sound belongs to the moment it ran out, not to the page load.
  all("button.timer").forEach(function (button) {
    var key = button.getAttribute("data-key");
    var total = parseFloat(button.getAttribute("data-seconds")) || 0;
    var left = total;
    var end = 0;
    var handle = null;
    var time = button.querySelector(".timer-time");
    var action = button.querySelector(".timer-action");

    function show(label) {
      time.textContent = clock(left);
      action.textContent = label;
    }

    function stop() {
      window.clearInterval(handle);
      handle = null;
      button.classList.remove("is-running");
    }

    function run() {
      handle = window.setInterval(tick, 250);
      button.classList.add("is-running");
      show("Pause");
    }

    function finish(sound) {
      stop();
      left = 0;
      button.classList.add("is-finished");
      show("Reset");
      if (sound) alarm();
    }

    function tick() {
      left = Math.max(0, (end - Date.now()) / 1000);
      if (left <= 0) finish(true);
      else show("Pause");
    }

    var kept = state[key];
    if (kept && typeof kept.ends === "number") {
      end = kept.ends;
      left = Math.max(0, (end - Date.now()) / 1000);
      if (left > 0) run();
      else finish(false);
    } else if (kept && typeof kept.left === "number") {
      left = kept.left;
      show("Resume");
    }

    button.addEventListener("click", function () {
      if (handle !== null) {
        stop();
        state[key] = { left: left };
        show("Resume");
      } else if (left <= 0) {
        left = total;
        button.classList.remove("is-finished");
        delete state[key];
        show("Start");
      } else {
        end = Date.now() + left * 1000;
        state[key] = { ends: end };
        run();
      }
      save();
    });
  });
})();
