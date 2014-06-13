v = dlmread('/home/user/code/pelops/karakum/build/comp.txt',' ');
plot(v(:,1),v(:,2),'o');
xlabel('Analytical (v)');
ylabel('BEM (v)');