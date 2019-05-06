[h1,hp] = boundedline(results.xvalues, results.rem.word, results.rem.wordSE,'b', results.xvalues, results.for.word, results.for.wordSE, 'c', 'transparency', 0.1)

figure

[h1,hp] = boundedline(results.xvalues, results.rem.face, results.rem.faceSE,'k', results.xvalues, results.for.face, results.for.faceSE, 'g','transparency', 0.1)
hold on
[h1,hp] = boundedline(results.xvalues, results.rem.scene, results.rem.sceneSE,'r', results.xvalues, results.for.scene, results.for.sceneSE, 'm','transparency', 0.1)
