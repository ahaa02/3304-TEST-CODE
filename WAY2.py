#모델 슬라이싱
mesh_min = aabb_bounds_a[4] #최솟값
mesh_max = aabb_bounds_a[5] #최댓값
mesh_len = mesh_max - mesh_min #모델 길이
mesh_nn = int(input("Enter the number of slices (2, 4, 8, 16, ...): "))

# 슬라이스 개수가 2의 n승인지 확인
if not (mesh_nn & (mesh_nn - 1) == 0 and mesh_nn != 0):
    raise ValueError("Number of slices must be a power of 2 (2, 4, 8, 16, ...).")
    
mesh_interval = mesh_len / float(mesh_nn) #나눠진 길이
mesh_mid = mesh_interval/2 #중간값

z_min_ls = [mesh_min + (mesh_interval*i) for i in range(0, mesh_nn)] #모든 최소 값들 리스트
z_max_ls = [mesh_min + (mesh_interval*i) for i in range(1, mesh_nn + 1)] #모든 최대 값들 리스트

bounds_df = pd.DataFrame(data=[], columns=['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']) #데이터 프레임

for i in range(mesh_nn): #모델 나누기
    bounds_df.loc[i] = aabb_bounds_a
bounds_df['zmin'] = z_min_ls
bounds_df['zmax'] = z_max_ls

for i in range(mesh_nn): #나눈 공간 박스화
    globals()['box_{}'.format(i)] = list(bounds_df.iloc[i])
    exec('extract_%s = pv_mesh.clip_box(box_%s, invert=False)' % (i, i))
    
# 2D 이미지 생성
image_data = []
output_dir = 'output_images'
os.makedirs(output_dir, exist_ok=True)

p = pv.Plotter()

for i in range(mesh_nn): #슬라이싱된 모델 출력
    color = (random.random(), random.random(), random.random())
    exec("n_cell = extract_%s.n_cells" % (i))
    evalCode = 'n_cell != 0'
    if eval(evalCode):
        #exec("pv.save_meshio('extract_%s.stl', extract_%s)" % (i, i))
        exec("p.add_mesh(extract_%s, color=%s)" % (i, color))
        exec("py_box_%s = pv.Box(bounds=box_%s, level=0, quads=True)" % (i, i))
        exec("p.add_mesh(py_box_%s, opacity=0.1, show_edges=True)" % (i))
        # 이미지 생성
        z_value = z_min_ls[i] + mesh_mid
        slice_2d = pv_mesh.slice(normal='z', origin=(0, 0, z_value))
        image_file = os.path.join(output_dir, f'slice_{i}.png')
        pv.plot(slice_2d, screenshot=image_file, off_screen=True)
        image_data.append(image_file)
        
P.show()

# 전개도의 정해진 길이 리스트 생성(2^n 순서)
def split_slices(n):
    slices = []
    power = 1
    slices.append(2)
    n -= 2
    while n > 0:
        power *= 2
        if n >= power:
            slices.append(power)
            n -= power
        else:
            slices.append(n)
            n = 0
    return slices
slice_groups = split_slices(mesh_nn)

# 슬라이스된 이미지 리스트 출력
print("Slice groups: ", slice_groups)
print("Image files: ", image_data)

# 전개도 배치 규칙
def arrange_images(slice_groups, image_data): 
    result = [None] * len(image_data)
    index = 0
    # 첫 번째 두 이미지는 result의 처음과 끝에 배치
    result[0] = image_data[index]
    result[-1] = image_data[index + 1]
    index += 2
    left_index = len(result) // 2 - 1
    right_index = len(result) // 2
    reverse = True
    for group_size in slice_groups[1:]:
        for i in range(group_size // 2):
            if reverse:
                result[left_index - i] = image_data[index + 1]
                result[right_index + i] = image_data[index]
            else:
                result[left_index - i] = image_data[index]
                result[right_index + i] = image_data[index + 1]
            index += 2
            reverse = not reverse
        left_index -= (group_size // 2)
        right_index += (group_size // 2)
    return result
# 배열된 이미지 결과
result = arrange_images(slice_groups, image_data)
print("Result: ", result)

from PIL import Image, ImageEnhance
import os

# 이미지들을 연결해서 하나의 이미지로 만드는 함수
def concatenate_images(image_paths, output_path):
    images = [Image.open(img_path).convert("RGBA") for img_path in image_paths if img_path]  # 이미지 경로가 유효한 경우에만 열기
    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)
    concatenated_image = Image.new('RGBA', (total_width, max_height), (0, 0, 0, 0))  # 투명 배경 이미지 생성
    x_offset = 0
    for img in images:
        concatenated_image.paste(img, (x_offset, 0), img)  # 이미지를 합성하며 투명한 부분은 그대로 유지
        x_offset += img.width
    return concatenated_image
# 검은 배경에 이미지를 놓는 함수
def overlay_on_black_background(image_paths, output_path):
    concatenated_image = concatenate_images(image_paths, output_path)
    # 검은 배경 이미지 생성
    black_background = Image.new('RGBA', concatenated_image.size, (0, 0, 0, 255))
    # 검은 배경 위에 이미지 붙이기
    final_image = Image.alpha_composite(black_background, concatenated_image)
    # 저장
    final_image.save(output_path)
    print(f"Final image saved at {output_path}")
    
# 배열된 이미지 경로 설정
output_folder = 'output_images'
os.makedirs(output_folder, exist_ok=True)
result = ['output_images\\slice_0.png', 'output_images\\slice_7.png', 'output_images\\slice_5.png',
          'output_images\\slice_2.png', 'output_images\\slice_3.png', 'output_images\\slice_4.png',
          'output_images\\slice_6.png', 'output_images\\slice_1.png']
