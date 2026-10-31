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

  // A thing carrying step keys — a section of this page, a page of the run — says which steps it
  // stands for and how far the bench got in them. The marks are a store's, never a second count.
  function stepKeys(item) {
    return (item.getAttribute("data-steps") || "").split(" ").filter(Boolean);
  }

  function ticked(marks, keys) {
    var done = 0;
    keys.forEach(function (key) {
      if (marks["step-" + key] === true) done += 1;
    });
    return done;
  }

  var groups = all("nav details[data-steps]");

  function showSections() {
    groups.forEach(function (group) {
      var keys = stepKeys(group);
      var label = group.querySelector(".section-progress");
      if (label) label.textContent = ticked(state, keys) + " of " + keys.length + " done";
    });
  }

  function currentSection() {
    return groups.filter(function (group) {
      var keys = stepKeys(group);
      return ticked(state, keys) < keys.length;
    })[0];
  }

  function keyOf(box) {
    return box.getAttribute("data-key");
  }

  function isTicked(box) {
    return box.checked;
  }

  // Into the page's store; `save` writes the store to the browser.
  function record(box) {
    if (box.checked) state[keyOf(box)] = true;
    else delete state[keyOf(box)];
  }

  // A step's own mark is its instructions together, so each count reads what the bench did
  // however the reader ticks: every instruction ticked ticks the step, one cleared clears it, and
  // the step ticked or cleared by hand does the same to all of them. A step with no instructions
  // is ticked by hand.
  var stepMarks = stepBoxes.map(function (box) {
    var step = box.closest(".step");
    return {
      box: box,
      instructions: step ? all('.instructions input[type="checkbox"][data-key]', step) : []
    };
  });

  function refresh() {
    var done = 0;
    stepMarks.forEach(function (mark) {
      var box = mark.box;
      var step = box.closest(".step");
      if (step) step.classList.toggle("is-done", box.checked);
      box.indeterminate = !box.checked && mark.instructions.some(isTicked);
      if (box.checked) done += 1;
    });
    if (progress) progress.textContent = done + " of " + steps(stepBoxes.length) + " done";
    showSections();
  }

  function settle(mark, changed) {
    if (changed === mark.box) {
      mark.instructions.forEach(function (box) { box.checked = mark.box.checked; });
    } else {
      mark.box.checked = mark.instructions.every(isTicked);
    }
    [mark.box].concat(mark.instructions).forEach(record);
    save();
    refresh();
  }

  boxes.forEach(function (box) {
    box.checked = state[keyOf(box)] === true;
  });

  // The store can hold a step's mark without its instructions' or the reverse. Nothing ticked is
  // lost: a ticked step stands for its instructions, and instructions all ticked for their step.
  // It is saved once mended, so a run's index, which reads step marks alone, counts it from then.
  var mended = false;
  stepMarks.forEach(function (mark) {
    if (!mark.instructions.length) return;
    if (!mark.box.checked && !mark.instructions.every(isTicked)) return;
    [mark.box].concat(mark.instructions).forEach(function (box) {
      if (!box.checked) mended = true;
      box.checked = true;
      record(box);
    });
  });
  if (mended) save();

  stepMarks.forEach(function (mark) {
    [mark.box].concat(mark.instructions).forEach(function (box) {
      box.addEventListener("change", function () { settle(mark, box); });
    });
  });
  refresh();

  // The section the bench is in — the first still holding an unticked step — opens on arrival,
  // and the rest close. A section opened by hand after that stays open. Where this never runs,
  // the page keeps the first section open, which is where a bench with no marks is.
  var current = currentSection();
  if (current) {
    groups.forEach(function (group) { group.open = group === current; });
  }

  // A run's index: how far the bench got in each protocol, read from that page's own store.
  // Every file:// page shares one store, and each page has a key of its own, so the index
  // reads the marks without the protocol pages writing anything twice.
  all("[data-page-key]").forEach(function (item) {
    // Each step of that page by the key it is addressed under, which a reworded title keeps.
    var keys = stepKeys(item);
    var label = item.querySelector(".page-progress");
    var marks;
    try {
      marks = JSON.parse(
        window.localStorage.getItem("liulab-protocol:" + item.getAttribute("data-page-key")) || "{}"
      );
    } catch (error) {
      return;
    }
    if (!marks || !label || !keys.length) return;
    var done = ticked(marks, keys);
    label.textContent = done + " of " + steps(keys.length) + " done";
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

  // A figure opens to its record's map, which loads the first time it is opened. Nothing is
  // remembered: an opened map is not a mark.
  all("details.opener").forEach(function (opener) {
    opener.addEventListener("toggle", function () {
      var frame = opener.querySelector("iframe[data-src]");
      if (opener.open && frame && !frame.getAttribute("src")) {
        frame.setAttribute("src", frame.getAttribute("data-src"));
      }
    });
  });

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
    var live = calculators(figure, function () { update(); });

    function update() {
      var reactions = Math.max(1, Math.floor(Number(input.value) || 1));
      var scale = reactions * (1 + overage);
      if (live) live.draw(scale);
      else {
        all("[data-ul]", figure).forEach(function (cell) {
          cell.textContent = number(parseFloat(cell.getAttribute("data-ul")) * scale);
        });
      }
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

  // Calculators. Arithmetic here is only over numbers render.py wrote, never a constant or a
  // formula of the package's.
  //
  // Every calculator's input is a reader's number: a field over the number the protocol gives,
  // its `data-plan`. The page keeps the reader's only where it differs, offers the protocol's
  // back on the button beside it, and holds the field still while anything in its step is
  // ticked, as a timer's time does while it runs. A held field says so, in render.py's words,
  // beside it and on itself. `changed` hears each value it takes.
  //
  // A field with no `data-plan` is one the protocol has no number for, such as a plate nobody
  // has counted: it opens empty, its value reads null until the reader types, and the button
  // empties it again.
  function readerNumber(field, changed) {
    var key = field.getAttribute("data-key");
    var plan = field.hasAttribute("data-plan")
      ? parseFloat(field.getAttribute("data-plan"))
      : null;
    var shown = field.defaultValue;
    var box = field.closest(".calc");
    var back = box.querySelector(".calc-plan");
    // Only a field something is pipetted against carries the word, and only it is held still.
    var held = box.querySelector(".calc-held");
    var least = parseFloat(field.getAttribute("data-least"));
    var step = field.closest(".step");
    var marks = step ? all('input[type="checkbox"][data-key]', step) : [];
    var value = allowed(state[key]) ? state[key] : plan;

    function locked() {
      return held !== null && marks.some(isTicked);
    }

    // What a field takes: more than nothing, or from its own `data-least` up where none is a
    // reading of its own.
    function allowed(given) {
      if (typeof given !== "number") return false;
      return isNaN(least) ? given > 0 : given >= least;
    }

    // What the reader typed: undefined for anything the field will not take. The text the page
    // opened with reads as the protocol's own, though a number is printed rounded.
    function typed() {
      var text = field.value.trim();
      if (text === shown) return plan;
      // The page groups a long number, so the reader may type the groups back.
      text = text.replace(/,/g, "");
      if (!/^\d*\.?\d+$/.test(text)) return undefined;
      return allowed(Number(text)) ? Number(text) : undefined;
    }

    function show() {
      field.value = value === plan ? shown : String(value);
      field.readOnly = locked();
      field.title = locked() ? held.textContent : "";
      if (held) held.hidden = !locked();
      back.hidden = value === plan || locked();
    }

    function take(to) {
      value = to;
      if (value !== plan) state[key] = value;
      else delete state[key];
      save();
      changed();
    }

    field.addEventListener("input", function () {
      var to = typed();
      if (to === undefined || locked()) return;
      take(to);
      back.hidden = value === plan;
    });
    // A slip leaves the number as it was, and the field says so once the reader moves on.
    field.addEventListener("change", show);
    back.addEventListener("click", function () { take(plan); show(); });
    if (step) {
      step.addEventListener("change", function (event) {
        if (event.target.type === "checkbox") show();
      });
    }
    show();
    return {
      plan: plan,
      value: function () { return value; }
    };
  }

  // A row the bench measures gives the volume that carries its amount, and the row it names
  // makes up the difference.
  function calculators(figure, redraw) {
    var rows = all("tbody > tr", figure);
    var measured = rows.filter(function (row) { return row.classList.contains("measured"); });
    if (!measured.length) return null;
    var planned = rows.map(function (row) { return parseFloat(row.getAttribute("data-rxn-ul")); });

    var calcs = measured.map(function (row) {
      return {
        row: row,
        at: rows.indexOf(row),
        fill: Number(row.getAttribute("data-fill")),
        nanograms: parseFloat(row.getAttribute("data-ng")),
        least: parseFloat(row.getAttribute("data-least-ul")),
        input: readerNumber(row.querySelector(".calc-value"), redraw)
      };
    });

    // Each row's volume for one reaction: the protocol's, or the one the typed concentration
    // gives, with what a row gained taken off the row that makes it up.
    function volumes() {
      var ul = planned.slice();
      calcs.forEach(function (one) {
        var value = one.input.value();
        if (value === one.input.plan) return;
        var volume = one.nanograms / value;
        ul[one.fill] -= volume - planned[one.at];
        ul[one.at] = volume;
      });
      return ul;
    }

    function draw(scale) {
      var ul = volumes();
      var fired = {};
      calcs.forEach(function (one) {
        var problem = one.row.getAttribute("data-too-dilute");
        if (problem && ul[one.fill] < 0) fired[problem] = true;
        problem = one.row.getAttribute("data-too-concentrated");
        if (problem && ul[one.at] < one.least) fired[problem] = true;
      });
      var each = 0;
      var whole = 0;
      var spent = false;
      rows.forEach(function (row, i) {
        // A row making up the difference cannot go below nothing: the reaction outgrows its
        // volume instead, and the total says by how much.
        var volume = Math.max(0, ul[i]);
        var short = ul[i] < 0;
        row.classList.toggle("is-over", short);
        spent = spent || short;
        row.querySelector(".one").textContent = number(volume);
        var mix = row.querySelector("[data-ul]");
        if (mix) {
          mix.textContent = number(volume * scale);
          each += volume;
        }
        whole += volume;
      });
      figure.classList.toggle("is-over", spent);
      var foot = figure.querySelector("tfoot");
      if (foot) {
        foot.querySelector(".one").textContent = number(whole);
        foot.querySelector("[data-ul]").textContent = number(each * scale);
      }
      all(".dispense-mix", figure).forEach(function (span) { span.textContent = number(each); });
      all(".dispense [data-row]", figure).forEach(function (span) {
        span.textContent = number(Math.max(0, ul[Number(span.getAttribute("data-row"))]));
      });
      all(".calc-warning", figure).forEach(function (warning) {
        warning.hidden = !fired[warning.getAttribute("data-problem")];
      });
    }

    return { draw: draw };
  }

  // A count the bench took: the control comes off it, the dilution scales it back up, and the
  // net is read against the floor render.py wrote. The floor itself is never worked out here,
  // and neither is anything the plan says about it. Until the plate is counted there is no net
  // and no verdict: the line keeps the dash render.py wrote, so nothing reads as a pass.
  all(".net-count").forEach(function (block) {
    var floor = parseFloat(block.getAttribute("data-floor"));
    var net = block.querySelector(".net-sum");
    var shown = block.querySelector(".net-value");
    var said = block.querySelector(".net-said");
    var verdict = block.querySelector(".net-verdict");
    var unread = shown.textContent;

    function part(name) {
      var input = block.querySelector("." + name + " .calc-value");
      return input ? readerNumber(input, function () { draw(); }) : null;
    }

    var counted = part("net-counted");
    var control = part("net-control");
    var dilution = part("net-dilution");

    function draw() {
      var plate = counted.value();
      var read = plate !== null;
      var total = read ? (plate - (control ? control.value() : 0)) * dilution.value() : 0;
      var clears = read && total >= floor;
      // A count is a number of things, so however many there are it is written out in full.
      shown.textContent = read
        ? total.toLocaleString("en-US", { maximumFractionDigits: 2 })
        : unread;
      said.hidden = !read;
      verdict.textContent = read ? net.getAttribute(clears ? "data-clears" : "data-short") : "";
      net.classList.toggle("is-short", read && !clears);
      all(".calc-warning", block).forEach(function (warning) {
        warning.hidden = !read || clears;
      });
    }

    draw();
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

  // What a reader types for a time: h:mm:ss, m:ss, or a bare number of minutes. Null for
  // anything else, or for no time at all, so a slip leaves the timer as it was.
  function parseTime(text) {
    var parts = String(text).trim().split(":");
    if (parts.length > 3) return null;
    var total = 0;
    for (var i = 0; i < parts.length; i += 1) {
      if (!/^\d+(\.\d+)?$/.test(parts[i].trim())) return null;
      total = total * 60 + Number(parts[i]);
    }
    if (parts.length === 1) total *= 60;
    return total > 0 ? total : null;
  }

  // A browser lets a page sound only once the reader has touched it, so the alarm plays on a
  // context a touch made or woke; one made when the time runs out can stay silent.
  var audio = null;

  function unlock() {
    try {
      var Context = window.AudioContext || window.webkitAudioContext;
      if (!audio && Context) audio = new Context();
      if (audio && audio.state === "suspended") audio.resume();
    } catch (error) {
      audio = null;
    }
  }

  function beep() {
    [0, 0.35, 0.7].forEach(function (offset) {
      var tone = audio.createOscillator();
      var gain = audio.createGain();
      tone.frequency.value = 880;
      gain.gain.value = 0.15;
      tone.connect(gain);
      gain.connect(audio.destination);
      tone.start(audio.currentTime + offset);
      tone.stop(audio.currentTime + offset + 0.2);
    });
  }

  // A context the browser holds back may be let go only at the reader's next touch, which is
  // no time for an alarm, so one let go more than two seconds late stays quiet.
  function alarm() {
    try {
      unlock();
      var asked = Date.now();
      if (audio && audio.state === "running") beep();
      else if (audio) {
        Promise.resolve(audio.resume()).then(function () {
          if (audio.state === "running" && Date.now() - asked < 2000) beep();
        }, function () {});
      }
    } catch (error) {
      /* no sound available */
    }
    if (navigator.vibrate) navigator.vibrate([200, 100, 200]);
  }

  // A tab out of sight is not heard by everyone, so the tab's own title says which ran out
  // until that timer is reset. A timer is named by its page's key and its own, since the panel
  // holds the run's other pages' timers too.
  var title = document.title;
  var here = document.body.getAttribute("data-protocol") || "";
  var ranOut = {};

  function timerId(page, key) {
    return page + "\n" + key;
  }

  function markTitle(id, label) {
    if (label) ranOut[id] = label;
    else delete ranOut[id];
    var labels = Object.keys(ranOut).map(function (one) { return ranOut[one]; });
    document.title = labels.length ? "Time up: " + labels.join(", ") + " · " + title : title;
  }

  // Pages of a run are often open in several tabs, and a timer that runs out sounds once among
  // them: each tab claims its deadline in the shared store and reads the claim back a moment
  // later, so two claiming at once settle on the last to write. A tab the browser lets sound
  // claims at once; one still waiting for a touch claims only if no such tab has. Every tab
  // flags its title.
  var tab = Math.random().toString(36).slice(2);
  var CLAIM_SETTLE_MS = 150;
  // The longest wait setTimeout honours; a later deadline is redrawn by the interval instead.
  var LONGEST_TIMEOUT_MS = 2147483647;
  // Redraw just past a deadline, so the entry reads done, and keep the page clear of the panel.
  var PAST_DEADLINE_MS = 50;
  var PANEL_GAP_PX = 24;
  var MUTED_TAB_WAIT_MS = 1000;
  var CLAIM_PREFIX = "liulab-protocol-alarm:";
  var DAY_MS = 86400000;

  function claimKey(page, key) {
    return CLAIM_PREFIX + page + ":" + key;
  }

  // A claim outlives its alarm so a late tab does not sound it again; after a day it is dropped.
  try {
    Object.keys(window.localStorage).forEach(function (name) {
      if (name.indexOf(CLAIM_PREFIX) !== 0) return;
      var ends = parseFloat(window.localStorage.getItem(name) || "");
      if (!(ends > Date.now() - DAY_MS)) window.localStorage.removeItem(name);
    });
  } catch (error) {
    /* nothing to tidy */
  }

  function claimed(claim, ends) {
    var held = window.localStorage.getItem(claim);
    return Boolean(held) && held.split(" ")[0] === String(ends);
  }

  function ring(page, key, label, ends) {
    markTitle(timerId(page, key), label);
    var claim = claimKey(page, key);
    var mine = ends + " " + tab;
    var audible = Boolean(audio) && audio.state === "running";
    function stake() {
      try {
        if (claimed(claim, ends)) return;
        window.localStorage.setItem(claim, mine);
      } catch (error) {
        alarm();
        return;
      }
      window.setTimeout(function () {
        var held = null;
        try {
          held = window.localStorage.getItem(claim);
        } catch (error) {
          /* unread: stay quiet rather than sound twice */
        }
        if (held === mine) alarm();
      }, CLAIM_SETTLE_MS);
    }
    if (audible) stake();
    else window.setTimeout(stake, MUTED_TAB_WAIT_MS);
  }

  function release(page, key) {
    markTitle(timerId(page, key), null);
    try {
      window.localStorage.removeItem(claimKey(page, key));
    } catch (error) {
      /* nothing was claimed */
    }
  }

  // A running timer keeps its deadline and a paused one the seconds it has left, so turning the
  // page or closing the tab does not lose an incubation. A deadline already past comes back
  // finished and silent: the sound belongs to the moment it ran out, not to the page load.
  // A time the reader typed is kept beside either, and is the timer's until they take the
  // plan's back; the plan's own is `data-seconds`.
  var timers = all(".timer[data-key]");
  // This page's timers as the panel reads them: each one's key, label, element and controls.
  var own = timers.map(function (timer) {
    var key = timer.getAttribute("data-key");
    var plan = parseFloat(timer.getAttribute("data-seconds")) || 0;
    var label = (timer.querySelector(".timer-label") || timer).textContent;
    var field = timer.querySelector(".timer-time");
    var action = timer.querySelector(".timer-action");
    var back = timer.querySelector(".timer-plan");
    var total = plan;
    var left = total;
    var end = 0;
    var handle = null;
    var due = null;

    function show(word) {
      field.value = clock(left);
      field.readOnly = handle !== null;
      action.textContent = word;
      if (back) back.hidden = total === plan || handle !== null;
    }

    // What the page keeps for this timer: the reader's time where it is not the plan's, and the
    // deadline or the seconds left where either is set.
    function keep(extra) {
      var held = extra || {};
      if (total !== plan) held.seconds = total;
      if (Object.keys(held).length) state[key] = held;
      else delete state[key];
      save();
    }

    function stop() {
      window.clearInterval(handle);
      window.clearTimeout(due);
      handle = null;
      due = null;
      timer.classList.remove("is-running");
    }

    // The interval redraws, and a hidden tab may run it once a minute; the deadline is a timeout
    // of its own, which a browser holds to the second. One past what a timeout can wait, about
    // 24 days, is left to the interval.
    function run() {
      var wait = Math.max(0, end - Date.now());
      handle = window.setInterval(tick, 250);
      if (wait < LONGEST_TIMEOUT_MS) due = window.setTimeout(function () { finish(true); }, wait);
      timer.classList.add("is-running");
      show("Pause");
    }

    function finish(sounded) {
      if (timer.classList.contains("is-finished")) return;
      stop();
      left = 0;
      timer.classList.add("is-finished");
      show("Reset");
      if (sounded) ring(here, key, label, end);
    }

    function tick() {
      left = Math.max(0, (end - Date.now()) / 1000);
      if (left <= 0) finish(true);
      else show("Pause");
    }

    // Back to a full count of `to`, standing ready.
    function ready(to) {
      stop();
      total = to;
      left = to;
      timer.classList.remove("is-finished");
      release(here, key);
      keep();
      show("Start");
    }

    // As the store holds it: on arrival, and again when another tab changes it.
    function restore() {
      var kept = state[key] || {};
      stop();
      timer.classList.remove("is-finished");
      total = typeof kept.seconds === "number" ? kept.seconds : plan;
      left = total;
      if (typeof kept.ends === "number") {
        end = kept.ends;
        left = Math.max(0, (end - Date.now()) / 1000);
        if (left > 0) run();
        else finish(false);
      } else if (typeof kept.left === "number") {
        left = kept.left;
        show("Resume");
      } else {
        show("Start");
      }
    }

    restore();

    action.addEventListener("click", function () {
      if (handle !== null) {
        left = Math.max(0, (end - Date.now()) / 1000);
        stop();
        keep({ left: left });
        show("Resume");
      } else if (left <= 0) {
        ready(total);
      } else {
        unlock();
        end = Date.now() + left * 1000;
        keep({ ends: end });
        run();
      }
    });

    // A typed time is the whole count, except on a paused timer, where it is what is left and
    // the count grows or shrinks by the difference: 60 minutes become 90 with the bench mid-way.
    field.addEventListener("change", function () {
      var typed = parseTime(field.value);
      var paused = !timer.classList.contains("is-finished") && left !== total;
      if (typed === null || handle !== null) {
        show(action.textContent);
      } else if (paused) {
        total = total - left + typed;
        left = typed;
        keep({ left: left });
        show("Resume");
      } else {
        ready(typed);
      }
    });

    if (back) back.addEventListener("click", function () { ready(plan); });

    return {
      key: key,
      label: label,
      element: timer,
      restore: restore,
      reset: function () { ready(total); }
    };
  });

  // The timer panel: every timer running, paused or finished, on this page and, on a page of a
  // run, on the run's other protocol pages, each read from its own page's store. The one to
  // finish soonest is on top; a finished one stays until it is dismissed or reset.
  var panel = document.querySelector(".timer-panel");
  var others = [];
  try {
    var listed = document.getElementById("run-timers");
    if (listed) others = JSON.parse(listed.textContent) || [];
    others = others.filter(function (page) {
      return page.key !== here && page.timers && page.timers.length;
    });
  } catch (error) {
    others = [];
  }

  function storeOf(page) {
    if (page === here) return state;
    try {
      return JSON.parse(window.localStorage.getItem("liulab-protocol:" + page) || "{}") || {};
    } catch (error) {
      return {};
    }
  }

  var RANK = { done: 0, running: 1, paused: 2 };

  function entries() {
    var now = Date.now();
    var list = [];

    function add(kept, entry) {
      if (!kept || typeof kept !== "object") return;
      if (typeof kept.ends === "number") {
        entry.ends = kept.ends;
        entry.left = Math.max(0, (kept.ends - now) / 1000);
        entry.status = entry.left > 0 ? "running" : "done";
      } else if (typeof kept.left === "number") {
        entry.left = kept.left;
        entry.status = "paused";
      } else {
        return;
      }
      entry.id = timerId(entry.page, entry.key);
      list.push(entry);
    }

    own.forEach(function (one) {
      add(state[one.key], {
        page: here, key: one.key, label: one.label, title: "",
        href: "#" + one.element.id, element: one.element, reset: one.reset
      });
    });
    others.forEach(function (page) {
      var store = storeOf(page.key);
      (page.timers || []).forEach(function (one) {
        add(store[one.key], {
          page: page.key, key: one.key, label: one.label, title: page.title, href: one.href
        });
      });
    });
    list.sort(function (a, b) { return RANK[a.status] - RANK[b.status] || a.left - b.left; });
    return list;
  }

  // Brought into view and lit for a moment, from the panel or from another page's link.
  function flag(element) {
    var still = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    element.scrollIntoView({ behavior: still ? "auto" : "smooth", block: "center" });
    element.classList.add("is-flagged");
    window.setTimeout(function () { element.classList.remove("is-flagged"); }, 2000);
  }

  // A timer on another page is set back to ready in that page's store, keeping the reader's
  // time; its page, if open, hears the change.
  function dismiss(entry) {
    if (entry.reset) {
      entry.reset();
    } else {
      var store = storeOf(entry.page);
      var kept = store[entry.key] || {};
      if (typeof kept.seconds === "number") store[entry.key] = { seconds: kept.seconds };
      else delete store[entry.key];
      try {
        window.localStorage.setItem("liulab-protocol:" + entry.page, JSON.stringify(store));
      } catch (error) {
        /* the store is out of reach, so the timer stays */
      }
      release(entry.page, entry.key);
    }
    draw();
  }

  function span(name, text) {
    var one = document.createElement("span");
    one.className = name;
    one.textContent = text;
    return one;
  }

  function rowFor(entry) {
    var row = document.createElement("li");
    var link = document.createElement("a");
    link.href = entry.href;
    link.appendChild(span("panel-label", entry.label));
    if (entry.title) link.appendChild(span("panel-page", entry.title));
    var button = document.createElement("button");
    button.type = "button";
    button.className = "panel-dismiss";
    button.textContent = "Dismiss";
    row.appendChild(link);
    row.appendChild(span("panel-time", ""));
    row.appendChild(button);
    if (entry.element) {
      link.addEventListener("click", function (event) {
        event.preventDefault();
        flag(entry.element);
      });
    }
    button.addEventListener("click", function () { dismiss(row.entry); });
    return row;
  }

  var rows = {};
  var seen = {};
  var soon = null;

  function draw() {
    var list = entries();
    var holder = panel.querySelector("ol");
    var now = {};
    list.forEach(function (entry) {
      now[entry.id] = entry.status;
      var row = rows[entry.id] || (rows[entry.id] = rowFor(entry));
      row.entry = entry;
      row.className = "is-" + entry.status;
      row.querySelector(".panel-time").textContent =
        entry.status === "done" ? "done" : clock(entry.left);
      row.querySelector(".panel-dismiss").hidden = entry.status !== "done";
      holder.appendChild(row);
      // This page's own timers sound from their steps; another page's sounds here once this tab
      // has seen it run out.
      if (entry.page !== here && entry.status === "done" && seen[entry.id] === "running") {
        ring(entry.page, entry.key, entry.label, entry.ends);
      }
    });
    Object.keys(rows).forEach(function (id) {
      if (now[id]) return;
      holder.removeChild(rows[id]);
      delete rows[id];
    });
    Object.keys(ranOut).forEach(function (id) {
      if (now[id] !== "done") markTitle(id, null);
    });
    seen = now;
    panel.hidden = !list.length;
    // The page's end scrolls clear of the panel, so the panel never covers the credit.
    document.body.style.paddingBottom = panel.hidden ? "" : panel.offsetHeight + PANEL_GAP_PX + "px";
    // Another page's next deadline is waited for to the second, as this page's own are.
    window.clearTimeout(soon);
    var next = list.filter(function (entry) {
      return entry.page !== here && entry.status === "running";
    })[0];
    if (next && next.left * 1000 < LONGEST_TIMEOUT_MS) {
      soon = window.setTimeout(draw, next.left * 1000 + PAST_DEADLINE_MS);
    }
  }

  if (panel && (own.length || others.length)) {
    draw();
    window.setInterval(draw, 500);
    // Another tab changed a store: this page's own timers take what it holds now.
    window.addEventListener("storage", function (event) {
      if (event.key === storageKey) {
        var fresh = load();
        own.forEach(function (one) {
          if (JSON.stringify(fresh[one.key]) === JSON.stringify(state[one.key])) return;
          if (fresh[one.key] === undefined) delete state[one.key];
          else state[one.key] = fresh[one.key];
          one.restore();
        });
      }
      draw();
    });
  }

  // Opened at a timer, as the panel on another page of the run links one.
  try {
    var target = document.getElementById(decodeURIComponent(window.location.hash.slice(1)));
    if (target && target.classList.contains("timer")) flag(target);
  } catch (error) {
    /* an address naming no timer */
  }

  // A timer still running after a reload sounds only once the page has been touched again.
  // A touch lets a page sound when it lifts, not when it lands.
  ["pointerup", "keydown"].forEach(function (kind) {
    document.addEventListener(kind, function () {
      if (document.querySelector(".timer.is-running, .timer-panel .is-running")) unlock();
    }, true);
  });
})();
