x = 1750:125:3000;
TWF1 = [0 0 0 34 70 93 100 100 100 100 100]/100;
TanhWF1 = [0 22 73 93 100 100 100 100 100 100 100]/100;
TWF2 = [0 0 0 1 9 17 44 69 78 91 93]/100;
TanhWF2 = [0 3 11 45 78 92 96 94 98 100 100]/100;
figure
plot(x, TWF1, 'b-o', x, TanhWF1, 'b-^', x, TWF2,'r--o', x, TanhWF2,'r--^')
set(gca,'XTick',x);
legend('TWF (Tanh Init)', 'TanhWF (Tanh Init)', 'TWF (Truncated Init)', ...
'TanhWF (Truncated Init)', 'location', 'southeast');
legend('boxoff');
xlabel('Number of measurements');
ylabel('Empirical success rate');