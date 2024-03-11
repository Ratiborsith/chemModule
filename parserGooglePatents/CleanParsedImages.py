import os
import shutil


def delete_every_second_image(folder_path):
    files = os.listdir(folder_path)
    image_files = [file for file in files if file.endswith(('.png', '.jpg', '.jpeg', '.gif', '.bmp'))]

    for i, file in enumerate(image_files):
        if i % 2 == 1:
            file_path = os.path.join(folder_path, file)
            os.remove(file_path)
            print(f"Deleted: {file_path}")


# Укажите путь к вашей папке с изображениями
folder_path = 'C:\\pycharmpro\\chemModule\\parserGooglePatents\\Patent_images\\US20230136730A1'
delete_every_second_image(folder_path)
