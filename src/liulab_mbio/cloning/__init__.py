"""The cloning methods: how an experiment joins its fragments into one plasmid.

One method is one package here, and `plan` is what every method's plan shares. Golden Gate and
Gibson assembly are the methods this package holds today. Nothing is re-exported: import a
method by module path, as `liulab_mbio.cloning.goldengate`.

`docs/adr/0007-cloning-methods.md` records why a method is a directory, here and in the skills,
and why the overhang rules sit below every method rather than inside one.
"""
