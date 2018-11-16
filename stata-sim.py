user_query = [x for x in input().split()]

if user_query[0] == 'sum':
    nums = [float(i) for i in user_query[1:]]
    user_result = sum(nums)
    print(user_result)
