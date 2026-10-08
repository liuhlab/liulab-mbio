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
    if (progress) progress.textContent = done + " of " + stepBoxes.length + " steps done";
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

  var clear = document.querySelector("button.clear");
  if (clear) {
    clear.addEventListener("click", function () {
      boxes.forEach(function (box) {
        box.checked = false;
        delete state[box.getAttribute("data-key")];
      });
      save();
      refresh();
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

    if (typeof state[key] === "number") input.value = String(state[key]);
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

  all("button.timer").forEach(function (button) {
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

    function tick() {
      left = Math.max(0, (end - Date.now()) / 1000);
      if (left <= 0) {
        stop();
        left = 0;
        button.classList.add("is-finished");
        show("Reset");
        alarm();
      } else {
        show("Pause");
      }
    }

    button.addEventListener("click", function () {
      if (handle !== null) {
        stop();
        show("Resume");
      } else if (left <= 0) {
        left = total;
        button.classList.remove("is-finished");
        show("Start");
      } else {
        end = Date.now() + left * 1000;
        handle = window.setInterval(tick, 250);
        button.classList.add("is-running");
        show("Pause");
      }
    });
  });
})();
