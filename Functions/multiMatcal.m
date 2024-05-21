function[ageprob] = multiMatcal(age, error, date_is)
resError = 200;
ageprob = zeros(55001, length(date_is));
    for i = date_is'
        [~,~,holder,~] = matcal(age(date_is(i))*1000, error(date_is(i))*1000,  'Marine20', 'CalBP','reserr', resError, 'plot', 0);
        ageprob(:,i) = holder(:,2);
    end
end