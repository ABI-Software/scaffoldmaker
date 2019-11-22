'''
Utility function for tests.
'''

def assertAlmostEqualList(testcase, actualList, expectedList, delta):
    assert len(actualList) == len(expectedList)
    for actual, expected in zip(actualList, expectedList):
        testcase.assertAlmostEqual(actual, expected, delta=delta)